classdef CoagulationRHS < handle
    % COAGULATIONRHS
    % ODE right-hand side for coagulation model.
    % Works for:
    %   - 0-D slab (use_column = false)
    %   - 1-D column (use_column = true)
    %
    % Uses:
    %   betas.b1, betas.b25   (b25 = b2 - b3 - b4 - b5)
    %   linear                (growth + sinking + other linear terms)
    %   optional linear disagg (0-D legacy)
    %   optional NPP source
    %
    % Nonlinear disaggregation and attenuation are applied
    % outside this RHS in CoagulationSimulation (run loop).

    properties
        betas          % BetaMatrices object
        linear         % Linear operator matrix
        disaggMinus    % Linear disaggregation (legacy)
        disaggPlus     % Linear disaggregation (legacy)
        config         % SimulationConfig
    end

    methods
        % --------------------------------------------------------------
        function obj = CoagulationRHS(betas, linear, disaggMinus, disaggPlus, config)
            obj.betas       = betas;
            obj.linear      = linear;
            obj.disaggMinus = disaggMinus;
            obj.disaggPlus  = disaggPlus;
            obj.config      = config;
        end

        % --------------------------------------------------------------
        function updateKernel(obj, newBetas)
            obj.betas = newBetas;
        end

        % --------------------------------------------------------------
        function dvdt = evaluate(obj, t, v)
            % EVALUATE  Compute dv/dt for coagulation + linear terms + NPP

            % 1) Make sure v is a column and non-negative
            v_col = v(:);
            v_pos = max(v_col, 0);

            % 2) Get coagulation scaling from epsilon, etc.
            [coag_scale, ~] = obj.getScalingFactors(t);

            % 3) Infer Ns from betas and Nz from v length
            Ns = obj.betas.getNumSections();   % number of size bins from kernels
            nY = numel(v_pos);

            if obj.config.use_column
                % ---- Column mode: v is [Nz * Ns x 1] ----
                Nz = nY / Ns;
                if abs(Nz - round(Nz)) > 1e-8
                    error('CoagulationRHS:StateSizeMismatch', ...
                        'Length(v) = %d not equal to Ns * Nz with Ns = %d', ...
                        nY, Ns);
                end
                Nz = round(Nz);

                % Arrange as [Nz x Ns]: rows = depth, cols = size class
                V_grid      = reshape(v_pos, Ns, Nz).';   % [Nz x Ns]
                dCoag_grid  = zeros(Nz, Ns);

                b25 = obj.betas.b25;
                b1  = obj.betas.b1;

                for k = 1:Nz
                    row       = V_grid(k, :);               % 1 x Ns
                    rowShift  = [0, row(1:end-1)];          % shifted PSD
                    term1     = row .* (row * b25);         % gain/loss combo
                    term2     = (row * b1) .* rowShift;     % super-diagonal
                    dCoag_grid(k, :) = term1 + term2;
                end

                dCoag_grid = dCoag_grid * coag_scale;      % scale by alpha(eps)
                dvdt_coag  = dCoag_grid.';                 % [Ns x Nz]
                dvdt_coag  = dvdt_coag(:);                 % back to column

                % Linear processes (sinking + growth + any other linear ops)
                dvdt_linear = obj.linear * v_pos;

                dvdt = dvdt_coag + dvdt_linear;

            else
                % ---- 0-D mode: v is [Ns x 1] ----
                if nY ~= Ns
                    error('CoagulationRHS:StateSizeMismatch0D', ...
                        '0-D state length = %d but Ns = %d', nY, Ns);
                end

                v_row    = v_pos.';                        % 1 x Ns
                v_shift  = [0, v_row(1:end-1)];           % shifted PSD
                term1    = v_row .* (v_row * obj.betas.b25);
                term2    = (v_row * obj.betas.b1) .* v_shift;
                dvdt_coag = (term1 + term2).' * coag_scale;

                dvdt = dvdt_coag + (obj.linear * v_pos);

                % Legacy linear disaggregation in 0-D, if nonlinear is OFF
                if ~obj.config.disagg_use_nonlinear ...
                        && ~isempty(obj.disaggMinus) && ~isempty(obj.disaggPlus)
                    dvdt = dvdt - (obj.disaggMinus * v_pos) + (obj.disaggPlus * v_pos);
                end
            end

            % 4) NPP Source term
            if obj.config.use_NPP
                source_vec = obj.getNPPSource(t, Ns, nY);
                dvdt = dvdt + source_vec;
            end

            % 5) Safety: NaN / Inf guard
            bad = ~isfinite(dvdt);
            if any(bad)
                dvdt(bad) = 0;
            end
        end

        % --------------------------------------------------------------
        function J = jacobian(obj, ~, ~)
            % For now, just return the linear operator as a Jacobian approx
            J = obj.linear;
        end

        function J_loc = computeLocalJacobian(obj, ~, ~)
            % Placeholder for local Jacobian; not used.
            J_loc = [];
        end

        % --------------------------------------------------------------
        function [coag_scale, eps_here] = getScalingFactors(obj, t)
            % Get epsilon(t) and corresponding coagulation scaling.

            if strcmpi(obj.config.epsilon_profile, 'observed') ...
                    && ~isempty(obj.config.epsilon_series)

                if ~isempty(obj.config.epsilon_time)
                    eps_here = interp1( ...
                        obj.config.epsilon_time(:), ...
                        obj.config.epsilon_series(:), ...
                        t, 'linear', 'extrap');
                else
                    eps_here = obj.config.epsilon;
                end
            else
                eps_here = obj.config.epsilon;
            end

            % Avoid zero epsilon
            eps_here = max(eps_here, 1e-14);

            % alpha(eps) scaling
            ratio     = eps_here / obj.config.epsilon_ref;
            alpha_eff = obj.config.alpha_base * (ratio ^ -obj.config.p_alpha);

            % clip alpha to reasonable range
            alpha_eff = min(max(alpha_eff, obj.config.alpha_clip_min), ...
                            obj.config.alpha_clip_max);

            coag_scale = alpha_eff;
        end

        % --------------------------------------------------------------
        function s = getNPPSource(obj, t, Ns, nY) %#ok<INUSD>
            % Build NPP source vector (same length as dvdt).
            %
            % For now:
            %   - constant rate NPP_rate
            %   - only in section 1
            %   - only in top layer if column mode

            rate = obj.config.NPP_rate;

            if rate == 0
                s = zeros(nY, 1);
                return;
            end

            if obj.config.use_column
                % Column: allocate full [Nz*Ns x 1] vector
                Nz = nY / Ns;
                if abs(Nz - round(Nz)) > 1e-8
                    error('CoagulationRHS:NPPSizeMismatch', ...
                        'NPP: length mismatch with Ns=%d, nY=%d', Ns, nY);
                end
                Nz = round(Nz);

                S_grid = zeros(Nz, Ns);
                % Put NPP in top layer, first size class
                S_grid(1, 1) = rate;

                s = S_grid.';     % [Ns x Nz]
                s = s(:);         % [Ns*Nz x 1]
            else
                % 0-D: only Ns entries
                s = zeros(Ns, 1);
                s(1) = rate;
            end
        end

        function validate(obj) %#ok<MANU>
            % Hook for future checks
        end

        function [term1, term2, term3, term4, term5] = rateTerms(obj, v, t) %#ok<INUSD>
            n = length(v);
            term1 = zeros(n,1);
            term2 = term1;
            term3 = term1;
            term4 = term1;
            term5 = term1;
        end
    end
end