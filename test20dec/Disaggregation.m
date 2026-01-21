classdef Disaggregation
    %DISAGGREGATION Turbulence-driven breakup operators (MASS-conserving)
    %
    % Legacy (kept for backward compatibility):
    %   Ynew = Disaggregation.apply(Y, grid, epsilon)
    %   Ynew = Disaggregation.apply(Y, grid, epsilon, params)
    %
    % New (used by CoagulationSimulation.run):
    %   Ynew = Disaggregation.applyWithScaling(Y, grid, eps_here, eps_ref, n_exp)
    %   Ynew = Disaggregation.applyWithScaling(..., params)
    %
    % IMPORTANT (NEW-2025-12-22):
    % This file now enforces MASS (particle volume) conservation:
    %   sum( Y .* Vbin ) is conserved, NOT sum(Y).
    %
    % Reason:
    %   Your export flux uses Vbin * Y, so conserving number breaks mass balance
    %   and can cause Relative Export Flux to explode (>1) when disagg is ON.
    %
    % 'params' for apply (legacy):
    %   .kappa (3.5)  – cutoff steepness in logistic
    %   .C0    (2e-3) – baseline size scale (cm)
    %   .B     (0.45) – ε exponent in r_max = C0 * ε^(-B)
    %   .p     (2.5)  – weight exponent for redistribution (r^-p)
    %
    % 'params' for applyWithScaling (overrides cfg when provided):
    %   .kmax_a (0.60)    – fraction of bins active at ε_ref
    %   .beta   (0.35)    – ε-sensitivity of active-tail width
    %   .fedge  (1/3)     – fraction of broken mass to immediate smaller bin (j-1)
    %   .pdist  (0)       – bias exponent for spreading remainder across [1..j-1]
    %
    % If grid.config exists, the following SimulationConfig fields are used:
    %   disagg_kmax_a, disagg_beta, disagg_frac_to_edge, disagg_redistribute_p

    methods (Static)

        % ============================
        % Legacy operator (same API, now MASS-conserving)
        % ============================
        function Ynew = apply(Yin, grid, epsilon, params)
            % Yin: state (vector), grid must supply radii
            % epsilon: W kg^-1

            if nargin < 4 || isempty(params)
                params = struct('kappa',3.5,'C0',2e-3,'B',0.45,'p',2.5);
            else
                if ~isfield(params,'kappa'), params.kappa = 3.5; end
                if ~isfield(params,'C0'),    params.C0    = 2e-3; end
                if ~isfield(params,'B'),     params.B     = 0.45; end
                if ~isfield(params,'p'),     params.p     = 2.5;  end
            end

            % columnize & guards
            Y0 = Yin(:);

            % ==========================================================
            % NEW-2025-12-22: compute bin "mass" (volume) weights Vbin
            % Use conserved radii if available (matches your flux code)
            % ==========================================================
            Vbin = Disaggregation.getBinVolumes_cm3(grid);

            % OLD (kept): number mass
            % mass0 = sum(Y0);

            % NEW: volume/mass in bin units
            mass0 = sum(Y0 .* Vbin);

            if ~isfinite(mass0) || mass0 <= 0 || ~isfinite(epsilon) || epsilon <= 0
                Ynew = Y0;  return;
            end

            % radii (cm) for the logistic threshold
            if ismethod(grid,'getFractalRadii')
                r = grid.getFractalRadii();  r = r(:);
            else
                error('Disaggregation:MissingRadii','grid.getFractalRadii() not found.');
            end

            % logistic mask that suppresses large bins as ε grows
            r_max = params.C0 * max(epsilon, realmin)^(-params.B);
            ratio = r ./ max(r_max, 1e-12);
            frag_factor = 1 ./ (1 + exp(params.kappa*(ratio - 1))); % ~1 below threshold, ~0 above
            Y_frag = Y0 .* frag_factor;

            % ==========================================================
            % OLD (kept): lost number
            % lost_mass = mass0 - sum(Y_frag);
            % ==========================================================

            % NEW: lost MASS (volume)
            lost_mass = mass0 - sum(Y_frag .* Vbin);
            if lost_mass < 0, lost_mass = 0; end

            Y = Y_frag;

            if lost_mass > 0
                % weights (use r^-p, but redistribute MASS, not number)
                w = r.^(-params.p);
                sw = sum(w);
                if sw > 0
                    add_mass = lost_mass * (w / sw);         % mass to add per bin
                    Y = Y + add_mass ./ max(Vbin, realmin);  % convert mass -> number
                else
                    % fallback: dump all lost MASS to smallest bin
                    Y(1) = Y(1) + lost_mass / max(Vbin(1), realmin);
                end
            end

            % ==========================================================
            % OLD (kept): exact number conservation
            % mass_new = sum(Y);
            % if mass_new > 0 && isfinite(mass_new)
            %     Y = Y * (mass0 / mass_new);
            % end
            % ==========================================================

            % NEW: exact MASS conservation
            mass_new = sum(Y .* Vbin);
            if mass_new > 0 && isfinite(mass_new)
                Y = Y * (mass0 / mass_new);
            end

            % hygiene
            Y(~isfinite(Y)) = 0;
            Y = max(Y, 0);

            Ynew = Y;
        end

        % ===========================================
        % New operator: ε-scaled, tail-focused breakup
        % NOW MASS-conserving in (Y .* Vbin)
        % ===========================================
        function Yout = applyWithScaling(Yin, grid, eps_here, eps_ref, n_exp, params)
            % eps_here – current ε(t) [W kg^-1]
            % eps_ref  – reference ε (e.g., 1e-6) for scaling
            % n_exp    – nonlinearity (0.6–0.9 gentle; 1.0 stronger)
            % params   – (optional) struct; see class header

            if nargin < 5 || isempty(n_exp),  n_exp  = 0.45; end
            if nargin < 4 || isempty(eps_ref), eps_ref = max(eps_here, realmin); end
            if nargin < 6, params = struct(); end

            % columnize & guards
            Yout = Yin(:);
            nsec = numel(Yout);
            if nsec < 3 || ~isfinite(eps_here) || eps_here <= 0
                return;
            end

            % ==========================================================
            % NEW-2025-12-22: Vbin for MASS conservation
            % ==========================================================
            Vbin = Disaggregation.getBinVolumes_cm3(grid);

            tiny = 1e-30;

            % pull defaults from cfg when available
            cfg = [];
            if isprop(grid,'config'); cfg = grid.config; end

            if isfield(params,'kmax_a')
                kmax_a = params.kmax_a;
            elseif ~isempty(cfg) && isprop(cfg,'disagg_kmax_a')
                kmax_a = cfg.disagg_kmax_a;
            else
                kmax_a = 0.60;
            end

            if isfield(params,'beta')
                beta = params.beta;
            elseif ~isempty(cfg) && isprop(cfg,'disagg_beta')
                beta = cfg.disagg_beta;
            else
                beta = 0.35;
            end

            if isfield(params,'fedge')
                fedge = params.fedge;
            elseif ~isempty(cfg) && isprop(cfg,'disagg_frac_to_edge')
                fedge = cfg.disagg_frac_to_edge;
            else
                fedge = 1/3;
            end
            fedge = min(max(fedge,0),1);

            if isfield(params,'pdist')
                pdist = params.pdist;
            elseif ~isempty(cfg) && isprop(cfg,'disagg_redistribute_p')
                pdist = cfg.disagg_redistribute_p;
            else
                pdist = 0;
            end

            % intensity of breakup relative to ε_ref (bounded)
            eps_rel   = max(eps_here, tiny) / max(eps_ref, tiny);
            intensity = 1 - eps_rel.^(-max(n_exp,0));      % →0 when ε≪ε_ref; ↑ as ε≫ε_ref
            if ~isfinite(intensity) || intensity <= 0
                return; % negligible breakup
            end
            intensity = min(max(intensity, 0), 0.95);

            % how many top (large) bins participate
            tail_scale = (1 ./ max(eps_rel, tiny)).^max(beta,0);
            kmax = max(1, round( min(nsec-1, kmax_a * tail_scale * nsec ) ));
            j1 = nsec - kmax + 1;
            tail_idx = j1:nsec;

            % gradation across the tail: last bin breaks most
            rloc = (1:numel(tail_idx)) / numel(tail_idx);
            local = 0.5 + 0.5 * rloc.^0.5;
            break_frac = intensity * local;
            break_frac = min(break_frac, 0.95);

            % ==========================================================
            % NEW-2025-12-22: store MASS before
            % ==========================================================
            mass_in = sum(Yin(:) .* Vbin);

            % redistribute tail mass (MASS-based transfers)
            for u = 1:numel(tail_idx)
                j = tail_idx(u);

                Nj = Yout(j);
                if Nj <= 0, continue; end

                % OLD (kept): break number
                % m_break = break_frac(u) * Nj;

                % NEW: break MASS in bin j
                mj_mass   = Nj * Vbin(j);
                m_break   = break_frac(u) * mj_mass;
                if m_break <= 0, continue; end

                % remove from source bin (convert mass -> number)
                Yout(j) = Nj - (m_break / max(Vbin(j), realmin));

                % 1) send part of broken MASS to immediate smaller bin
                if j > 1
                    Yout(j-1) = Yout(j-1) + (fedge * m_break) / max(Vbin(j-1), realmin);
                else
                    % safety
                    Yout(j) = Yout(j) + (fedge * m_break) / max(Vbin(j), realmin);
                end

                % 2) spread remaining MASS across smaller bins [1..j-1]
                rem = (1 - fedge) * m_break;
                if rem > 0 && j > 1
                    k = (1:(j-1))';
                    if pdist > 0
                        w = (j - k) .^ pdist;
                    else
                        w = ones(j-1,1);
                    end
                    sw = sum(w);
                    if sw > 0
                        add_mass = rem * (w / sw); % mass to each smaller bin
                        Yout(k) = Yout(k) + add_mass ./ max(Vbin(k), realmin);
                    else
                        Yout(j-1) = Yout(j-1) + rem / max(Vbin(j-1), realmin);
                    end
                end
            end

            % hygiene
            Yout(~isfinite(Yout)) = 0;
            Yout = max(Yout, 0);

            % ==========================================================
            % NEW-2025-12-22: exact MASS conservation fix
            % ==========================================================
            mass_out = sum(Yout .* Vbin);
            dm = mass_in - mass_out;   % mass error

            if abs(dm) > 0 && isfinite(dm)
                kfix = find(Yout>0, 1, 'first');
                if isempty(kfix), kfix = 1; end
                Yout(kfix) = Yout(kfix) + dm / max(Vbin(kfix), realmin);
                Yout = max(Yout, 0);
            end
        end

        % ==============================================================
        % Helper: bin volumes in cm^3 (constant factors cancel in ratios)
        % ==============================================================        
        function Vbin = getBinVolumes_cm3(grid)
            % Prefer conserved radii (matches your export flux code)
            if ismethod(grid,'getConservedRadii')
                rv = grid.getConservedRadii(); rv = rv(:);
            elseif ismethod(grid,'getFractalRadii')
                rv = grid.getFractalRadii();  rv = rv(:);
            else
                error('Disaggregation:MissingRadii','Need grid.getConservedRadii() or getFractalRadii().');
            end
            Vbin = (4/3) * pi * max(rv,0).^3;     % cm^3 (if rv is cm)
            Vbin(~isfinite(Vbin) | Vbin<=0) = realmin;
        end

    end
end