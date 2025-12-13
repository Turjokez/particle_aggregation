classdef CoagulationRHS < handle
    % COAGULATIONRHS ODE right-hand side for coagulation equations
    %
    % coagulation + growth + sinking + disaggregation.
    %


    properties
        betas;          % BetaMatrices object, contains coagulation parameters
        linear;         % Linear matrix (growth - sinking) OR column operator
        disaggMinus;    % Disaggregation loss matrix (Ns x Ns)
        disaggPlus;     % Disaggregation gain matrix (Ns x Ns)
        config;         % SimulationConfig object, stores configuration parameters

        % ==============================================================
        % NEW-2025-12-11: epsilon(t) storage for turbulence forcing
        % ==============================================================
        eps_time = [];     % NEW-2025-12-11: time vector [d]
        eps_vals = [];     % NEW-2025-12-11: epsilon(t) [e.g., m^2 s^{-3}]
        eps_const = [];    % NEW-2025-12-11: fallback constant epsilon
    end

    methods
        function obj = CoagulationRHS(betas, linear, disaggMinus, disaggPlus, config)
            % COAGULATIONRHS Constructor
            obj.betas       = betas;
            obj.linear      = linear;
            obj.disaggMinus = disaggMinus;
            obj.disaggPlus  = disaggPlus;
            obj.config      = config;

            % NEW-2025-12-11: store default constant epsilon from config (if present)
            if isprop(config, 'epsilon_const')
                obj.eps_const = config.epsilon_const;
            else
                obj.eps_const = [];
            end
        end

        % NEW-allow CoagulationSimulation to provide eps(t)
        % ==============================================================
        function setEpsilonSeries(obj, t, eps_vals)
            obj.eps_time = t(:);
            obj.eps_vals = eps_vals(:);
        end

         
        % NEW-compatibility wrapper for MATLAB ODE solvers
        % Many wrappers expect a method named rhs(t,y). Keep evaluate() too.

        function dvdt = rhs(obj, t, v)
            dvdt = obj.evaluate(t, v);
        end

        function dvdt = evaluate(obj, t, v)
            % EVALUATE Evaluate ODE right-hand side
            %
            % 0-D slab:
            %   v = [Ns x 1]
            %
            % 1-D column (NEW-2025-12-11):
            %   v = [Ns*Nz x 1], internally reshaped to N = [Ns x Nz]

            Ns = obj.config.n_sections;

            % NEW-2025-12-11: get epsilon(t) for this time (for info/use)

            if ~isempty(obj.eps_time) && ~isempty(obj.eps_vals)
                eps_now = interp1(obj.eps_time, obj.eps_vals, t, 'linear', 'extrap'); %#ok<NASGU>
            elseif ~isempty(obj.eps_const)
                eps_now = obj.eps_const; %#ok<NASGU>
            else
                eps_now = []; %#ok<NASGU>
            end

            % NEW- optional control for negative handling
            %
            % IMPORTANT FIX (closure-safe):
            %   - Linear operator MUST see RAW state v (v_lin)
            %   - Coag/disagg can see clamped state (v_nl)

            clip_negative = true;  % NEW-2025-12-12 (default)
            if isprop(obj.config, 'clip_negative')
                clip_negative = logical(obj.config.clip_negative);
            end

            % ==========================================================
            % NEW-split state for closure safety
            % ==========================================================
            v_lin = v;                       % RAW (for linear operator)
            if clip_negative
                v_nl = max(v, 0);            % NONNEG (for coag/disagg)
            else
                v_nl = v;
            end

            % ==========================================================
            % NEW-detect growth_mode for Adrian PP tests
            % ==========================================================
            growth_mode = "shift";
            if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                growth_mode = lower(string(obj.config.growth_mode));
            end
            mu = obj.config.growth;   % Adrian mu stored here

            % ==========================================================
            % NEW- prefer mass-conserving beta form if available
            % If BetaAssembler created b25 to be conservative, use it.
            % ==========================================================
            has_b25 = false;
            try
                has_b25 = isprop(obj.betas,'b25') && ~isempty(obj.betas.b25);
            catch
                has_b25 = false;
            end

            % ==========================================================
            % NEW-Column branch
            % ==========================================================
            if isprop(obj.config, 'use_column') && obj.config.use_column
                Nz = obj.config.getNumLayers();

                % Safety check
                if numel(v) ~= Ns * Nz
                    error('Column RHS: state length mismatch. Expected Ns*Nz = %d, got %d.', Ns*Nz, numel(v));
                end

                % Reshape NONLINEAR state into Ns x Nz
                N = reshape(v_nl, [Ns, Nz]);   % Ns rows (size), Nz columns (depth)

                % Preallocate coag+disagg tendency in Ns x Nz
                dN = zeros(Ns, Nz);

                % --------- apply coag + disaggregation layer-by-layer ---------
                c3 = obj.config.c3; %#ok<NASGU>
                c4 = obj.config.c4; %#ok<NASGU>

                for k = 1:Nz
                    vk = N(:, k);  % Ns x 1

                    % safe vector for coag/disagg only
                    vk_safe = max(vk, 0);

                    vk_r     = vk_safe';               % 1 x Ns
                    vk_shift = [0, vk_r(1:Ns-1)];      % 1 x Ns

                    % ======================================================
                    % NEW-2025-12-13: CONSERVATION-SAFE coag form
                    % Prefer b25 if present, otherwise fallback to b2-b3-b4-b5
                    % ======================================================
                    % (re-check has_b25 safely; cheap and robust)
                    has_b25_local = false;
                    try
                        has_b25_local = isprop(obj.betas,'b25') && ~isempty(obj.betas.b25);
                    catch
                        has_b25_local = has_b25;
                    end

                    if has_b25_local
                        term1 = vk_r * obj.betas.b25;      % 1 x Ns
                        term1 = vk_r .* term1;             % 1 x Ns
                    else
                        term_gain2 = vk_r * obj.betas.b2;  % 1 x Ns
                        term_gain2 = vk_r .* term_gain2;

                        term_loss3 = vk_r * obj.betas.b3;
                        term_loss3 = vk_r .* term_loss3;

                        term_loss4 = vk_r * obj.betas.b4;
                        term_loss4 = vk_r .* term_loss4;

                        term_loss5 = vk_r * obj.betas.b5;
                        term_loss5 = vk_r .* term_loss5;

                        term1 = (term_gain2 - term_loss3 - term_loss4 - term_loss5); % 1xNs
                    end

                    term2 = vk_r * obj.betas.b1;           % 1 x Ns
                    term2 = term2 .* vk_shift;             % 1 x Ns

                    dvk = (term1 + term2)';                % Ns x 1

                    % Disaggregation (matrix form)
                    dvk = dvk - obj.disaggMinus * vk_safe + obj.disaggPlus * vk_safe;

                    dN(:, k) = dvk;
                end

                % --------- add linear operator (includes vertical transfer + surface growth) ---------
                % CRITICAL FIX:
                %   linear operator MUST see RAW state v_lin (closure-safe)
                term3_flat = obj.linear * v_lin;              % (Ns*Nz) x 1
                term3 = reshape(term3_flat, [Ns, Nz]);        % Ns x Nz

                dN = dN + term3;

                % ==========================================================
                % NEW-2025-12-14: Adrian PP term (surface layer only)
                % ==========================================================
                if growth_mode == "pp"
                    % state ordering: v = [N(:,1); N(:,2); ...]
                    idx_N1_surface = 1;  % N(1,1)
                    dN(1,1) = dN(1,1) + mu * v_lin(idx_N1_surface);
                end

                dvdt = dN(:);
                return;
            end

            % ==========================================================
            % 0-D slab branch
            % ==========================================================
            n_sections = length(v); %#ok<NASGU>  % Number of sections (legacy)

            % Nonlinear uses v_nl, linear uses v_lin
            v_r     = v_nl';
            v_shift = [0, v_r(1:Ns-1)];

            % ======================================================
            % NEW-2025-12-13: CONSERVATION-SAFE coag form (0-D)
            % ======================================================
            if has_b25
                term1 = v_r * obj.betas.b25;
                term1 = v_r .* term1;
            else
                term_gain2 = v_r * obj.betas.b2;
                term_gain2 = v_r .* term_gain2;

                term_loss3 = v_r * obj.betas.b3;
                term_loss3 = v_r .* term_loss3;

                term_loss4 = v_r * obj.betas.b4;
                term_loss4 = v_r .* term_loss4;

                term_loss5 = v_r * obj.betas.b5;
                term_loss5 = v_r .* term_loss5;

                term1 = (term_gain2 - term_loss3 - term_loss4 - term_loss5);
            end

            term2 = v_r * obj.betas.b1;
            term2 = term2 .* v_shift;

            % CRITICAL FIX:
            % linear operator MUST see RAW state v_lin (closure-safe)
            term3 = obj.linear * v_lin;

            % Disaggregation (use nonnegative nonlinear state)
            term4 = -obj.disaggMinus * max(v_nl,0);
            term5 =  obj.disaggPlus  * max(v_nl,0);

            dvdt = (term1 + term2)' + term3 + term4 + term5;

            % ==========================================================
            % NEW-2025-12-14: Adrian PP term (0-D)
            % ==========================================================
            if growth_mode == "pp"
                dvdt(1) = dvdt(1) + mu * v_lin(1);
            end
        end

        function J = jacobian(obj, t, v) %#ok<INUSD>
            % JACOBIAN Evaluate analytical Jacobian
            %
            % IMPORTANT NEW-2025-12-12:
            % Default: SAFE approximation (linear + disagg). Full Jacobian optional.

            use_full_jacobian = false; % NEW-2025-12-12 default
            if isprop(obj.config,'use_full_jacobian') && ~isempty(obj.config.use_full_jacobian)
                use_full_jacobian = logical(obj.config.use_full_jacobian);
            end

            % ==========================================================
            % Column mode: safe Jacobian approx
            % ==========================================================
            if isprop(obj.config, 'use_column') && obj.config.use_column
                Ns = obj.config.n_sections;
                Nz = obj.config.getNumLayers();

                Iz = speye(Nz);
                D  = (-obj.disaggMinus + obj.disaggPlus);

                J = obj.linear + kron(Iz, D);

                % ----------------------------------------------------------
                % NEW-2025-12-14: add PP Jacobian term if growth_mode='pp'
                % ----------------------------------------------------------
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

            % ==========================================================
            % 0-D mode: safe Jacobian approx by default
            % ==========================================================
            if ~use_full_jacobian
                J = obj.linear - obj.disaggMinus + obj.disaggPlus;
                return;
            end

            % ==========================================================
            % OLD full Jacobian (kept, but wrapped in try/catch)
            % ==========================================================
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

                % NEW-2025-12-14: add PP Jacobian for 0-D if needed
                growth_mode = "shift";
                if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                    growth_mode = lower(string(obj.config.growth_mode));
                end
                if growth_mode == "pp"
                    mu = obj.config.growth;
                    J(1,1) = J(1,1) + mu;
                end

                return;
            end  % <<< ADDED: closes catch/try
        end      % <<< ADDED: closes jacobian() method

        function [term1, term2, term3, term4, term5] = rateTerms(obj, v)
            % RATETERMS Evaluate individual rate terms for diagnostics
            %
            % IMPORTANT FIX:
            %   term3 (linear) must use RAW v, not clamped v

            Ns = obj.config.n_sections;

            clip_negative = true;
            if isprop(obj.config, 'clip_negative')
                clip_negative = logical(obj.config.clip_negative);
            end

            v_lin = v;
            if clip_negative
                v_nl = max(v,0);
            else
                v_nl = v;
            end

            % NEW-2025-12-13: prefer b25 if present
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

                % CRITICAL FIX:
                term3 = obj.linear * v_lin;

                Dm = obj.disaggMinus;
                Dp = obj.disaggPlus;

                T4 = zeros(Ns, Nz);
                T5 = zeros(Ns, Nz);

                for k = 1:Nz
                    vk = N(:, k);
                    vk_safe = max(vk, 0);

                    T4(:, k) = -Dm * vk_safe;
                    T5(:, k) =  Dp * vk_safe;
                end

                term4 = T4(:);
                term5 = T5(:);
                return;
            end

            % ---------------- 0-D ----------------
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

            % CRITICAL FIX:
            term3 = obj.linear * v_lin;

            term4 = -obj.disaggMinus * max(v_nl,0);
            term5 =  obj.disaggPlus  * max(v_nl,0);
        end

        function validate(obj)
            % VALIDATE Validate RHS configuration
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