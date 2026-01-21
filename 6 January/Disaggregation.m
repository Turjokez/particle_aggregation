classdef Disaggregation
    %DISAGGREGATION Turbulence-driven breakup operators (BIOVOLUME-conserving)
    %
    % Legacy (kept for backward compatibility):
    %   Ynew = Disaggregation.apply(Y, grid, epsilon)
    %   Ynew = Disaggregation.apply(Y, grid, epsilon, params)
    %
    % New (used by CoagulationSimulation.run):
    %   Ynew = Disaggregation.applyWithScaling(Y, grid, eps_here, eps_ref, n_exp)
    %   Ynew = Disaggregation.applyWithScaling(..., params)
    %
    % IMPORTANT (FIX-2026-01-12):
    % Your model state Y is already in BIOVOLUME concentration per bin
    % (validated by: inventory from state = inventory from output).
    %
    % Therefore the disaggregation operator must conserve:
    %   sum(Y)   (per depth layer)
    % NOT sum(Y .* Vbin).
    %
    % If you conserve Y.*Vbin while Y is already biovolume, you will create
    % artificial mass (exactly what you observed: dv_disagg sums positive).
    %
    % 'params' for apply (legacy):
    %   .kappa (3.5)  – cutoff steepness in logistic
    %   .C0    (2e-3) – baseline size scale (cm)
    %   .B     (0.45) – ε exponent in r_max = C0 * ε^(-B)
    %   .p     (2.5)  – redistribution bias (r^-p) for WHERE mass goes
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
        % Legacy operator (same API, now BIOVOLUME-conserving)
        % ============================
        function Ynew = apply(Yin, grid, epsilon, params)
            % Yin: state vector (biovolume per bin)
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

            % BIOVOLUME in this layer / vector
            mass0 = sum(Y0);

            if ~isfinite(mass0) || mass0 <= 0 || ~isfinite(epsilon) || epsilon <= 0
                Ynew = Y0;
                return;
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

            % Apply mask (removes some biovolume from "too-large" bins)
            Y_frag = Y0 .* frag_factor;

            % Lost BIOVOLUME that must be redistributed to smaller bins
            lost_mass = mass0 - sum(Y_frag);
            if lost_mass < 0, lost_mass = 0; end

            Y = Y_frag;

            if lost_mass > 0
                % Redistribute lost BIOVOLUME across bins with a radius bias
                w = r.^(-params.p);
                sw = sum(w);
                if sw > 0
                    Y = Y + lost_mass * (w / sw);
                else
                    % fallback: dump all lost mass to smallest bin
                    Y(1) = Y(1) + lost_mass;
                end
            end

            % Exact BIOVOLUME conservation (renormalize)
            mass_new = sum(Y);
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
        % NOW BIOVOLUME-conserving in sum(Y)
        % ===========================================
        function Yout = applyWithScaling(Yin, grid, eps_here, eps_ref, n_exp, params)
            % eps_here – current ε(t) [W kg^-1]
            % eps_ref  – reference ε (e.g., 1e-6) for scaling
            % n_exp    – nonlinearity (0.6–0.9 gentle; 1.0 stronger)
            % params   – (optional) struct; see class header

            if nargin < 5 || isempty(n_exp),   n_exp   = 0.45; end
            if nargin < 4 || isempty(eps_ref), eps_ref = max(eps_here, realmin); end
            if nargin < 6, params = struct(); end

            % columnize & guards
            Yout = Yin(:);
            nsec = numel(Yout);
            if nsec < 3 || ~isfinite(eps_here) || eps_here <= 0
                return;
            end

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
            intensity = 1 - eps_rel.^(-max(n_exp,0));      % ->0 when eps<<eps_ref; increases when eps>>eps_ref
            if ~isfinite(intensity) || intensity <= 0
                return; % negligible breakup
            end
            intensity = min(max(intensity, 0), 0.95);

            % tail widens as eps increases
            tail_scale = (max(eps_rel, 1)).^max(beta,0);
            kmax = max(1, round( min(nsec-1, kmax_a * tail_scale * nsec ) ));

            j1 = nsec - kmax + 1;
            tail_idx = j1:nsec;

            % gradation across the tail: last bin breaks most
            rloc = (1:numel(tail_idx)) / numel(tail_idx);
            local = 0.5 + 0.5 * rloc.^0.5;
            break_frac = intensity * local;
            break_frac = min(break_frac, 0.95);

            % BIOVOLUME before (for exact conservation check/fix)
            mass_in = sum(Yin(:));

            % redistribute tail mass (BIOVOLUME-based transfers)
            for u = 1:numel(tail_idx)
                j = tail_idx(u);

                mj = Yout(j);
                if mj <= 0, continue; end

                % break BIOVOLUME from bin j
                m_break = break_frac(u) * mj;
                if m_break <= 0, continue; end

                % remove from source bin
                Yout(j) = mj - m_break;

                % 1) send part of broken mass to immediate smaller bin
                if j > 1
                    Yout(j-1) = Yout(j-1) + fedge * m_break;
                else
                    % safety
                    Yout(j) = Yout(j) + fedge * m_break;
                end

                % 2) spread remaining mass across smaller bins [1..j-1]
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
                        Yout(k) = Yout(k) + rem * (w / sw);
                    else
                        Yout(j-1) = Yout(j-1) + rem;
                    end
                end
            end

            % hygiene
            Yout(~isfinite(Yout)) = 0;
            Yout = max(Yout, 0);

            % Exact BIOVOLUME conservation fix (should be tiny)
            mass_out = sum(Yout);
            dm = mass_in - mass_out;

            if isfinite(dm) && abs(dm) > 0
                kfix = find(Yout>0, 1, 'first');
                if isempty(kfix), kfix = 1; end
                Yout(kfix) = max(Yout(kfix) + dm, 0);
            end
        end

    end
end