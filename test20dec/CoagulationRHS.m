classdef CoagulationRHS < handle
    % COAGULATIONRHS ODE right-hand side for coagulation equations
    %
    % coagulation + growth + sinking + disaggregation.
    %
    % NEW-2025-12-21 (CRITICAL UNITS FIX):
    % Solver time is in DAYS, so RHS must return dv/dt in [per day].
    % Any scaling by day_to_sec or z_max is WRONG here.
    %
    % Symptom of wrong scaling:
    %   ||A*v|| / ||rhs(t,v)|| ~ day_to_sec/z_max  (~1e3)
    %
    % Fix:
    %   keep any old scaling lines, but COMMENTED OUT.
    %
    % NEW-2025-12-21:
    % Optional debug printing (cfg.debug_rhs_units = true)

    properties
        betas;          % BetaMatrices object, contains coagulation parameters
        linear;         % Linear matrix (growth - sinking) OR column operator
        disaggMinus;    % Disaggregation loss matrix (Ns x Ns)
        disaggPlus;     % Disaggregation gain matrix (Ns x Ns)
        config;         % SimulationConfig object, stores configuration parameters

        % ==============================================================
        % NEW-2025-12-11: epsilon(t) storage for turbulence forcing
        % ==============================================================
        eps_time = [];     % time vector [d]
        eps_vals = [];     % epsilon(t) [vector OR Nt x Nz]
        eps_const = [];    % fallback constant epsilon

        % ==============================================================
        % NEW-2025-12-20: store grid so Disaggregation.applyWithScaling can use radii
        % ==============================================================
        grid = [];         % DerivedGrid (optional; required for new disagg)

        % ==============================================================
        % NEW-2025-12-21: cached z-grid for eps_fun / eps(t,z) (optional)
        % ==============================================================
        z_cache = [];      % cache z-centers from config (m)
    end

    methods
        function obj = CoagulationRHS(betas, linear, disaggMinus, disaggPlus, config, varargin)
            % COAGULATIONRHS Constructor
            obj.betas       = betas;
            obj.linear      = linear;
            obj.disaggMinus = disaggMinus;
            obj.disaggPlus  = disaggPlus;
            obj.config      = config;

            % NEW-2025-12-21: optional signature print (TEMP; keep commented)
            % disp('*** Using CoagulationRHS NEW-2025-12-21 (NO SCALING) ***');

            % NEW-2025-12-11: store default constant epsilon from config (if present)
            if isprop(config, 'epsilon_const')
                obj.eps_const = config.epsilon_const;
            else
                obj.eps_const = [];
            end

            % ==============================================================
            % NEW-2025-12-20: allow passing grid as 6th input (no breaking older calls)
            % Usage:
            %   CoagulationRHS(betas, linear, Dm, Dp, cfg, grid)
            % ==============================================================
            if nargin >= 6 && ~isempty(varargin)
                try
                    obj.grid = varargin{1};
                catch
                    obj.grid = [];
                end
            end

            % ==============================================================
            % NEW-2025-12-21: cache depth grid if possible (safe)
            % ==============================================================
            try
                if ismethod(config,'getZ')
                    obj.z_cache = config.getZ();
                end
            catch
                obj.z_cache = [];
            end
        end

        % ==============================================================
        % allow CoagulationSimulation to provide eps(t)
        % ==============================================================
        function setEpsilonSeries(obj, t, eps_vals)
            obj.eps_time = t(:);

            % OLD (kept):
            % obj.eps_vals = eps_vals(:);

            % NEW-2025-12-21: allow eps_vals to be vector OR Nt x Nz matrix
            if isvector(eps_vals)
                obj.eps_vals = eps_vals(:);

                % enforce matching length for vector case
                Nt = numel(obj.eps_time);
                Ne = numel(obj.eps_vals);
                n  = min(Nt, Ne);
                obj.eps_time = obj.eps_time(1:n);
                obj.eps_vals = obj.eps_vals(1:n);
            else
                % if eps matrix came as Nz x Nt, transpose it
                try
                    Nt = numel(obj.eps_time);
                    if size(eps_vals,2) == Nt && size(eps_vals,1) ~= Nt
                        eps_vals = eps_vals.';
                    end
                catch
                end
                obj.eps_vals = eps_vals; % Nt x Nz
            end
        end

        % ==============================================================
        % optional setter for grid (safe)
        % ==============================================================
        function setGrid(obj, gridObj)
            obj.grid = gridObj;
        end

        % compatibility wrapper for MATLAB ODE solvers
        function dvdt = rhs(obj, t, v)
            dvdt = obj.evaluate(t, v);
        end

        % ==============================================================
        % helper to get eps(t) from stored series (scalar eps)
        % ==============================================================
        function eps_now = getEpsScalar(obj, t)
            eps_now = [];
            if ~isempty(obj.eps_time) && ~isempty(obj.eps_vals) && isvector(obj.eps_vals)
                Nt = numel(obj.eps_time);
                Ne = numel(obj.eps_vals);
                n  = min(Nt, Ne);
                tt = obj.eps_time(1:n);
                ee = obj.eps_vals(1:n);
                try
                    eps_now = interp1(tt(:), ee(:), t, 'linear', 'extrap');
                catch
                    eps_now = [];
                end
            elseif ~isempty(obj.eps_const)
                eps_now = obj.eps_const;
            end

            % clamp nonnegative if numeric
            try
                if ~isempty(eps_now) && isnumeric(eps_now) && isfinite(eps_now)
                    eps_now = max(eps_now, 0);
                end
            catch
            end
        end

        function dvdt = evaluate(obj, t, v)
            % EVALUATE Evaluate ODE right-hand side
            %
            % 0-D slab:
            %   v = [Ns x 1]
            %
            % 1-D column:
            %   v = [Ns*Nz x 1], reshaped to N = [Ns x Nz]

            Ns = obj.config.n_sections;

            % epsilon(t) (optional)
            eps_now = obj.getEpsScalar(t); %#ok<NASGU>

            % negative handling
            clip_negative = true;
            if isprop(obj.config, 'clip_negative')
                clip_negative = logical(obj.config.clip_negative);
            end

            % ==========================================================
            % Disaggregation switch
            % ==========================================================
            use_new_disagg = false;
            if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                if logical(obj.config.enable_disagg)
                    if isprop(obj.config,'disagg_apply_in') && ~isempty(obj.config.disagg_apply_in)
                        use_new_disagg = strcmpi(string(obj.config.disagg_apply_in), "rhs");
                    else
                        use_new_disagg = true;
                    end
                end
            end
            if use_new_disagg && isempty(obj.grid)
                use_new_disagg = false;
            end

            % split state
            v_lin = v;
            if clip_negative
                v_nl = max(v, 0);
            else
                v_nl = v;
            end

            % growth mode
            growth_mode = "shift";
            if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                growth_mode = lower(string(obj.config.growth_mode));
            end
            mu = obj.config.growth;

            % prefer mass-conserving beta form if available
            has_b25 = false;
            try
                has_b25 = isprop(obj.betas,'b25') && ~isempty(obj.betas.b25);
            catch
                has_b25 = false;
            end

            % ==========================================================
            % Column branch
            % ==========================================================
            if isprop(obj.config, 'use_column') && obj.config.use_column
                Nz = obj.config.getNumLayers();

                if numel(v) ~= Ns * Nz
                    error('Column RHS: state length mismatch. Expected %d, got %d.', Ns*Nz, numel(v));
                end

                N  = reshape(v_nl, [Ns, Nz]);   % Ns x Nz
                dN = zeros(Ns, Nz);

                % z-grid
                zc = [];
                try
                    if ~isempty(obj.z_cache) && numel(obj.z_cache) == Nz
                        zc = obj.z_cache(:);
                    elseif ismethod(obj.config,'getZ')
                        zc = obj.config.getZ();
                        zc = zc(:);
                        obj.z_cache = zc;
                    end
                catch
                    zc = [];
                end
                if isempty(zc)
                    zc = ((1:Nz) - 0.5) * obj.config.dz;
                    zc = zc(:);
                end

                % epsilon input options
                has_eps_fun = false;
                try
                    if isprop(obj.config,'eps_fun') && ~isempty(obj.config.eps_fun)
                        has_eps_fun = isa(obj.config.eps_fun,'function_handle');
                    end
                catch
                    has_eps_fun = false;
                end

                has_eps_matrix = false;
                try
                    if isprop(obj.config,'epsilon_time') && isprop(obj.config,'epsilon_series')
                        if ~isempty(obj.config.epsilon_time) && ~isempty(obj.config.epsilon_series)
                            if ~isvector(obj.config.epsilon_series) && size(obj.config.epsilon_series,2) == Nz
                                has_eps_matrix = true;
                            end
                        end
                    end
                catch
                    has_eps_matrix = false;
                end

                for k = 1:Nz
                    vk = N(:, k);

                    eps_here = [];

                    if has_eps_fun
                        try
                            eps_here = obj.config.eps_fun(t, zc(k));
                        catch
                            eps_here = [];
                        end
                    end

                    if isempty(eps_here) && has_eps_matrix
                        try
                            tt = obj.config.epsilon_time(:);
                            Em = obj.config.epsilon_series; % Nt x Nz
                            eps_here = interp1(tt, Em(:,k), t, 'linear', 'extrap');
                        catch
                            eps_here = [];
                        end
                    end

                    if isempty(eps_here)
                        eps_here = obj.getEpsScalar(t);
                    end

                    % new disagg (optional)
                    if use_new_disagg
                        if ~isempty(eps_here) && isfinite(eps_here) && eps_here > 0
                            eps_ref = eps_here;
                            if isprop(obj.config,'eps_ref') && ~isempty(obj.config.eps_ref)
                                eps_ref = obj.config.eps_ref;
                            end
                            n_exp = 0.45;
                            if isprop(obj.config,'disagg_n_exp') && ~isempty(obj.config.disagg_n_exp)
                                n_exp = obj.config.disagg_n_exp;
                            end
                            vk = Disaggregation.applyWithScaling(vk, obj.grid, eps_here, eps_ref, n_exp);
                        end
                    end

                    vk_safe  = max(vk, 0);
                    vk_r     = vk_safe';
                    vk_shift = [0, vk_r(1:Ns-1)];

                    if has_b25
                        term1 = vk_r * obj.betas.b25;
                        term1 = vk_r .* term1;
                    else
                        term_gain2 = vk_r * obj.betas.b2;  term_gain2 = vk_r .* term_gain2;
                        term_loss3 = vk_r * obj.betas.b3;  term_loss3 = vk_r .* term_loss3;
                        term_loss4 = vk_r * obj.betas.b4;  term_loss4 = vk_r .* term_loss4;
                        term_loss5 = vk_r * obj.betas.b5;  term_loss5 = vk_r .* term_loss5;
                        term1 = (term_gain2 - term_loss3 - term_loss4 - term_loss5);
                    end

                    term2 = vk_r * obj.betas.b1;
                    term2 = term2 .* vk_shift;

                    dvk = (term1 + term2)';

                    % legacy matrix disagg
                    if ~use_new_disagg
                        dvk = dvk - obj.disaggMinus * vk_safe + obj.disaggPlus * vk_safe;
                    end

                    dN(:, k) = dvk;
                end

                % linear operator uses RAW state
                term3_flat = obj.linear * v_lin;

                % OLD (kept): WRONG scaling (DO NOT USE)
                % try
                %     term3_flat = term3_flat * (obj.config.z_max / obj.config.day_to_sec);
                % catch
                % end

                % NEW-2025-12-21: NO scaling here
                term3 = reshape(term3_flat, [Ns, Nz]);
                dN = dN + term3;

                % PP source (surface only)
                if growth_mode == "pp"
                    dN(1,1) = dN(1,1) + mu * v_lin(1);
                end

                dvdt = dN(:);

                % OLD (kept): WRONG global scaling (DO NOT USE)
                % dvdt = dvdt * (obj.config.z_max / obj.config.day_to_sec);

                % ==========================================================
                % NEW-2025-12-21: optional RHS units debug check
                % ==========================================================
                try
                    if isprop(obj.config,'debug_rhs_units') && logical(obj.config.debug_rhs_units)
                        dv_lin = obj.linear * v_lin;
                        dv_rhs = dvdt;

                        fprintf('[DEBUG RHS UNITS] t=%.4f d | ||A*v||_inf=%.3e | ||rhs||_inf=%.3e | ratio=%.3f\n', ...
                            t, norm(dv_lin, inf), norm(dv_rhs, inf), norm(dv_lin, inf)/max(norm(dv_rhs, inf), eps));
                    end
                catch
                end

                return;
            end

            % ==========================================================
            % 0-D slab branch
            % ==========================================================
            v_r     = v_nl';
            v_shift = [0, v_r(1:Ns-1)];

            if has_b25
                term1 = v_r * obj.betas.b25;
                term1 = v_r .* term1;
            else
                term_gain2 = v_r * obj.betas.b2;  term_gain2 = v_r .* term_gain2;
                term_loss3 = v_r * obj.betas.b3;  term_loss3 = v_r .* term_loss3;
                term_loss4 = v_r * obj.betas.b4;  term_loss4 = v_r .* term_loss4;
                term_loss5 = v_r * obj.betas.b5;  term_loss5 = v_r .* term_loss5;
                term1 = (term_gain2 - term_loss3 - term_loss4 - term_loss5);
            end

            term2 = v_r * obj.betas.b1;
            term2 = term2 .* v_shift;

            term3 = obj.linear * v_lin;

            % OLD (kept): WRONG scaling (DO NOT USE)
            % try
            %     term3 = term3 * (obj.config.z_max / obj.config.day_to_sec);
            % catch
            % end

            if ~use_new_disagg
                term4 = -obj.disaggMinus * max(v_nl,0);
                term5 =  obj.disaggPlus  * max(v_nl,0);
            else
                term4 = 0 * v_lin;
                term5 = 0 * v_lin;
            end

            dvdt = (term1 + term2)' + term3 + term4 + term5;

            if growth_mode == "pp"
                dvdt(1) = dvdt(1) + mu * v_lin(1);
            end

            % OLD (kept): WRONG global scaling (DO NOT USE)
            % dvdt = dvdt * (obj.config.z_max / obj.config.day_to_sec);
        end

        function J = jacobian(obj, t, v) %#ok<INUSD>
            use_full_jacobian = false;
            if isprop(obj.config,'use_full_jacobian') && ~isempty(obj.config.use_full_jacobian)
                use_full_jacobian = logical(obj.config.use_full_jacobian);
            end

            if isprop(obj.config, 'use_column') && obj.config.use_column
                Ns = obj.config.n_sections;
                Nz = obj.config.getNumLayers();

                Iz = speye(Nz);

                use_new_disagg = false;
                if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                    if logical(obj.config.enable_disagg)
                        if isprop(obj.config,'disagg_apply_in') && ~isempty(obj.config.disagg_apply_in)
                            use_new_disagg = strcmpi(string(obj.config.disagg_apply_in), "rhs");
                        else
                            use_new_disagg = true;
                        end
                    end
                end
                if use_new_disagg && isempty(obj.grid)
                    use_new_disagg = false;
                end

                if use_new_disagg
                    D = sparse(Ns, Ns);
                else
                    D = (-obj.disaggMinus + obj.disaggPlus);
                end

                J = obj.linear + kron(Iz, D);

                growth_mode = "shift";
                if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                    growth_mode = lower(string(obj.config.growth_mode));
                end
                if growth_mode == "pp"
                    mu = obj.config.growth;
                    J(1,1) = J(1,1) + mu;
                end
                return;
            end

            if ~use_full_jacobian
                J = obj.linear - obj.disaggMinus + obj.disaggPlus;
                return;
            end

            try
                n_sections = length(v);
                v_r = v';

                v_mat   = v_r(ones(1, n_sections), :);
                v_shift = [zeros(n_sections, 1), v_mat(:, 1:end-1)];

                bnet = (obj.betas.b2 - obj.betas.b3 - obj.betas.b4 - obj.betas.b5);
                term1_tmp = v_r * bnet;
                term1     = diag(term1_tmp) + diag(v_r) * bnet;

                term2a = v_r * obj.betas.b1;
                term2a = diag(term2a(2:end), -1);

                term2b = diag(obj.betas.b1, 1);
                term2b = term2b' .* v_r(1:end-1);
                term2b = diag(term2b, -1);

                term2c = diag(v_r(2:end), -1) .* bnet';
                term2  = term2a + term2b + term2c;

                term3a = obj.betas.b1  .* v_shift;
                term3b = bnet .* v_mat;

                term3  = (term3a + term3b)';
                term3  = triu(term3, 2) + tril(term3, -1);

                J = term1 + term2 + term3 + obj.linear - obj.disaggMinus + obj.disaggPlus;

            catch ME
                warning('Full Jacobian failed (%s). Falling back to safe Jacobian.', ME.message);
                J = obj.linear - obj.disaggMinus + obj.disaggPlus;

                growth_mode = "shift";
                if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                    growth_mode = lower(string(obj.config.growth_mode));
                end
                if growth_mode == "pp"
                    mu = obj.config.growth;
                    J(1,1) = J(1,1) + mu;
                end
                return;
            end
        end

        function [term1, term2, term3, term4, term5] = rateTerms(obj, v)
            Ns = obj.config.n_sections;

            clip_negative = true;
            if isprop(obj.config, 'clip_negative')
                clip_negative = logical(obj.config.clip_negative);
            end

            use_new_disagg = false;
            if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                if logical(obj.config.enable_disagg)
                    if isprop(obj.config,'disagg_apply_in') && ~isempty(obj.config.disagg_apply_in)
                        use_new_disagg = strcmpi(string(obj.config.disagg_apply_in), "rhs");
                    else
                        use_new_disagg = true;
                    end
                end
            end
            if use_new_disagg && isempty(obj.grid)
                use_new_disagg = false;
            end

            v_lin = v;
            if clip_negative
                v_nl = max(v,0);
            else
                v_nl = v;
            end

            has_b25 = false;
            try
                has_b25 = isprop(obj.betas,'b25') && ~isempty(obj.betas.b25);
            catch
                has_b25 = false;
            end

            if isprop(obj.config, 'use_column') && obj.config.use_column
                Nz = obj.config.getNumLayers();

                if numel(v) ~= Ns * Nz
                    error('rateTerms: Column state length mismatch. Expected %d, got %d.', Ns*Nz, numel(v));
                end

                N = reshape(v_nl, [Ns, Nz]);

                T1 = zeros(Ns, Nz);
                T2 = zeros(Ns, Nz);

                for k = 1:Nz
                    vk = N(:, k);
                    vk_safe = max(vk, 0);

                    vk_r     = vk_safe';
                    vk_shift = [0, vk_r(1:Ns-1)];

                    if has_b25
                        a1 = vk_r * obj.betas.b25;
                        a1 = vk_r .* a1;
                    else
                        bnet = (obj.betas.b2 - obj.betas.b3 - obj.betas.b4 - obj.betas.b5);
                        a1 = vk_r * bnet;
                        a1 = vk_r .* a1;
                    end

                    a2 = vk_r * obj.betas.b1;
                    a2 = a2 .* vk_shift;

                    T1(:, k) = a1';
                    T2(:, k) = a2';
                end

                term1 = T1(:);
                term2 = T2(:);

                term3 = obj.linear * v_lin;

                Dm = obj.disaggMinus;
                Dp = obj.disaggPlus;

                T4 = zeros(Ns, Nz);
                T5 = zeros(Ns, Nz);

                for k = 1:Nz
                    vk = N(:, k);
                    vk_safe = max(vk, 0);

                    if ~use_new_disagg
                        T4(:, k) = -Dm * vk_safe;
                        T5(:, k) =  Dp * vk_safe;
                    else
                        T4(:, k) = 0 * vk_safe;
                        T5(:, k) = 0 * vk_safe;
                    end
                end

                term4 = T4(:);
                term5 = T5(:);
                return;
            end

            % 0-D
            v_r     = v_nl';
            v_shift = [0, v_r(1:Ns-1)];

            if has_b25
                term1 = v_r * obj.betas.b25;
                term1 = v_r .* term1;
            else
                bnet = (obj.betas.b2 - obj.betas.b3 - obj.betas.b4 - obj.betas.b5);
                term1 = v_r * bnet;
                term1 = v_r .* term1;
            end

            term2 = v_r * obj.betas.b1;
            term2 = term2 .* v_shift;

            term3 = obj.linear * v_lin;

            if ~use_new_disagg
                term4 = -obj.disaggMinus * max(v_nl,0);
                term5 =  obj.disaggPlus  * max(v_nl,0);
            else
                term4 = 0 * v_lin;
                term5 = 0 * v_lin;
            end
        end

        function validate(obj)
            if isempty(obj.betas) || isempty(obj.linear)
                error('RHS not properly initialized');
            end

            Ns = obj.config.n_sections;

            if isprop(obj.config, 'use_column') && obj.config.use_column
                Nz = obj.config.getNumLayers();
                expected = Ns * Nz;
                if size(obj.linear, 1) ~= expected
                    error('Linear matrix dimension mismatch (column). Expected %d, got %d.', expected, size(obj.linear, 1));
                end
            else
                if size(obj.linear, 1) ~= Ns
                    error('Linear matrix dimension mismatch (0-D). Expected %d, got %d.', Ns, size(obj.linear, 1));
                end
            end

            if obj.betas.getNumSections() ~= Ns
                error('Beta matrices dimension mismatch');
            end
        end
    end
end