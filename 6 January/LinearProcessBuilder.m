% ==============================================================
% FILE: LinearProcessBuilder.m
% ==============================================================
classdef LinearProcessBuilder < handle
    %LINEARPROCESSBUILDER Builds linear operators for growth, sinking, and disaggregation
    %
    % (your full header kept as-is)
    %
    % NEW-2025-12-21 (SAFETY):
    % If using NEW disaggregation in RHS (cfg.enable_disagg + disagg_apply_in='rhs'),
    % then legacy disaggregation matrices should be ZERO to prevent any accidental
    % double counting elsewhere.
    %
    % NEW-2026-01-12:
    % - Add optional switches:
    %     cfg.enable_sinking (default true if missing)
    %     cfg.enable_linear  (default true if missing)
    % - If sinking disabled: return S = 0 matrix and sink_rate=0

    methods (Static)

        function G = growthMatrix(config, grid) %#ok<INUSD>
            n_sections = config.n_sections;

            growth_mode = 'shift';
            if isprop(config, 'growth_mode') && ~isempty(config.growth_mode)
                growth_mode = lower(string(config.growth_mode));
            end

            growth_loss = zeros(n_sections, 1);
            growth_gain = zeros(n_sections-1, 1);

            if config.gro_sec > 0
                growth_loss(config.gro_sec:n_sections-1) = -1;
                growth_gain(config.gro_sec:end) = 2;
            end

            G0_shift = diag(growth_loss) + diag(growth_gain, -1);

            % OLD (kept):
            % G0(1,1) = 1;
            % G0 = config.growth * G0;

            G0_shift(1,1) = 1;                 % legacy boundary trick (kept)
            G0_shift = config.growth * G0_shift;

            G0_pp = sparse(n_sections, n_sections);

            % OLD (kept but disabled):
            % G0_pp(1,1) = config.growth;

            % NEW-2025-12-14: PP handled in CoagulationRHS, so keep ZERO:
            G0_pp(1,1) = 0;

            if growth_mode == "pp"
                G0 = G0_pp;
            else
                G0 = G0_shift;
            end

            if ~isprop(config, 'use_column') || ~config.use_column
                G = sparse(G0);
                return;
            end

            Nz = config.getNumLayers();

            G = sparse(n_sections*Nz, n_sections*Nz);
            G(1:n_sections, 1:n_sections) = sparse(G0);
        end


        function [S, sink_rate, settling_vel_cmday] = sinkingMatrix(config, grid)

            % ==========================================================
            % NEW-2026-01-12: Optional sink disable switch
            % ==========================================================
            sink_on = true;
            try
                if isprop(config,'enable_sinking') && ~isempty(config.enable_sinking)
                    sink_on = logical(config.enable_sinking);
                end
            catch
                sink_on = true;
            end

            if ~sink_on
                Ns = config.n_sections;
                settling_vel_cmday = zeros(Ns,1);
                sink_rate = zeros(Ns,1);

                if ~isprop(config,'use_column') || ~config.use_column
                    S = sparse(Ns, Ns);
                    return;
                end

                Nz = config.getNumLayers();
                S = sparse(Ns*Nz, Ns*Nz);
                return;
            end

            % ------------------- original code (kept) -------------------
            fractal_radius   = grid.getFractalRadii();
            conserved_radius = grid.getConservedRadii();

            settling_vel = SettlingVelocityService.velocity(fractal_radius, conserved_radius, grid.setcon);
            settling_vel_cmday = settling_vel * config.day_to_sec;

            dz_eff = config.dz;          % [m]
            dz_eff_cm = dz_eff * 100;    % [cm]

            sink_rate = settling_vel_cmday / dz_eff_cm;   % [1/day]

            if ~isprop(config, 'use_column') || ~config.use_column
                S = spdiags(sink_rate, 0, numel(sink_rate), numel(sink_rate));
                return;
            end

            sinking_form = "flux";
            if isprop(config, 'sinking_form') && ~isempty(config.sinking_form)
                sinking_form = lower(string(config.sinking_form));
            end

            Nz = config.getNumLayers();
            Ns = config.n_sections;

            W = spdiags(sink_rate, 0, Ns, Ns);
            Ntot = Ns * Nz;

            if sinking_form == "loop"
                S = sparse(Ntot, Ntot);

                for k = 1:Nz
                    rows = (k-1)*Ns + (1:Ns);
                    cols = rows;

                    S(rows, cols) = S(rows, cols) + W;

                    if k >= 2
                        cols_above = (k-2)*Ns + (1:Ns);
                        S(rows, cols_above) = S(rows, cols_above) - W;
                    end
                end
            else
                Iz = speye(Nz);
                shift_down = spdiags(ones(Nz,1), -1, Nz, Nz);
                Tz = Iz - shift_down;

                S = kron(Tz, W);
            end

            if isprop(config, 'debug_sinking') && ~isempty(config.debug_sinking) && config.debug_sinking
                vtest = abs(rand(Ntot,1));
                dtest = -(S * vtest);

                Nmat = reshape(vtest, [Ns, Nz]);
                Nbottom = Nmat(:, end);
                export_loss_rate = sum(sink_rate .* Nbottom);

                colsum = sum(dtest);
                fprintf('[debug_sinking] sum(dN/dt) = %.6e,  -export_loss_rate = %.6e\n', colsum, -export_loss_rate);
            end
        end


        function [Dminus, Dplus] = disaggregationMatrices(config)
            n_sections = config.n_sections;

            % ==========================================================
            % NEW-2025-12-21 SAFETY:
            % If NEW disagg is applied in RHS, force legacy matrices to ZERO
            % ==========================================================
            try
                if isprop(config,'enable_disagg') && ~isempty(config.enable_disagg) && logical(config.enable_disagg)
                    if isprop(config,'disagg_apply_in') && ~isempty(config.disagg_apply_in)
                        if strcmpi(string(config.disagg_apply_in), "rhs")
                            Dminus = sparse(n_sections, n_sections);
                            Dplus  = sparse(n_sections, n_sections);
                            return;
                        end
                    end
                end
            catch
                % do nothing; fall through to legacy behavior
            end

            disagg_mode = "legacy";
            if isprop(config,'disagg_mode') && ~isempty(config.disagg_mode)
                disagg_mode = lower(string(config.disagg_mode));
            end

            Dminus = zeros(n_sections);
            Dplus  = zeros(n_sections);

            if disagg_mode == "legacy"
                % --------- BEGIN ORIGINAL LEGACY CODE (DO NOT CHANGE) ----------
                Dminus = zeros(n_sections);
                Dplus  = zeros(n_sections);
                if n_sections > 2
                    idx = 2:(n_sections-1);
                    for k = idx
                        Dminus(k, k)   = config.c3 * config.c4^k;
                        Dplus(k, k-1)  = config.c3 * config.c4^(k+1);
                    end
                end
                % --------- END ORIGINAL LEGACY CODE (DO NOT CHANGE) ------------
                return;
            end

            if n_sections > 2
                idx = 2:(n_sections-1);
                for k = idx
                    a = config.c3 * config.c4^k;
                    Dminus(k, k)   = a;
                    Dplus(k, k+1)  = a * config.c4;
                end
            end

            % Optional: use sparse (safe, but not required)
            % Dminus = sparse(Dminus);
            % Dplus  = sparse(Dplus);
        end


        function L = linearMatrix(config, grid)

            % ==========================================================
            % NEW-2026-01-12: Optional linear disable switch
            % If enable_linear=false => L=0
            % ==========================================================
            lin_on = true;
            try
                if isprop(config,'enable_linear') && ~isempty(config.enable_linear)
                    lin_on = logical(config.enable_linear);
                end
            catch
                lin_on = true;
            end

            if ~lin_on
                Ns = config.n_sections;
                if ~isprop(config,'use_column') || ~config.use_column
                    L = sparse(Ns, Ns);
                    return;
                end
                Nz = config.getNumLayers();
                L = sparse(Ns*Nz, Ns*Nz);
                return;
            end

            G = LinearProcessBuilder.growthMatrix(config, grid);
            [S, ~, ~] = LinearProcessBuilder.sinkingMatrix(config, grid);
            L = G - S;
        end

    end
end
