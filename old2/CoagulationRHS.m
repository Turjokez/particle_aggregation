classdef CoagulationRHS < handle
    % COAGULATIONRHS ODE right-hand side for coagulation equations
    
    properties
        betas;          % BetaMatrices object, contains coagulation parameters
        linear;         % Linear matrix (growth - sinking)
        disaggMinus;    % Disaggregation loss matrix
        disaggPlus;     % Disaggregation gain matrix
        config;         % SimulationConfig object, stores configuration parameters
    end

    methods
        function obj = CoagulationRHS(betas, linear, disaggMinus, disaggPlus, config)
            % COAGULATIONRHS Constructor
            % Initializes the CoagulationRHS object with the required parameters
            obj.betas = betas;        % Set beta matrices
            obj.linear = linear;      % Set linear matrix
            obj.disaggMinus = disaggMinus;  % Set disaggregation loss matrix
            obj.disaggPlus = disaggPlus;    % Set disaggregation gain matrix
            obj.config = config;      % Set configuration object
        end

        % dynamic kernel update support
        function updateKernel(obj, newBetas)
            % UPDATEKERNEL Replace the coagulation kernels used by the RHS
            % newBetas: BetaMatrices-like struct updated each timestep
            obj.betas = newBetas;
        end

        function dvdt = evaluate(obj, t, v)
            % EVALUATE Evaluate ODE right-hand side
            % t = time (unused but required by ODE solver)
            % v = concentration vector
            % Returns: dvdt = rate of change (the right-hand side of ODE)

            n_sections = length(v);  % Number of sections (length of concentration vector)

            % Match legacy behavior: clamp negative values to eps only for calculations
            v_pos = max(v, eps);  % Ensure non-negative concentrations

            v_r = v_pos';  % Convert the concentration vector to a row vector for matrix operations
            v_shift = [0, v_r(1:n_sections-1)];  % Shift vector for previous section alignment

            % --- Coagulation terms ---
            term1 = v_r * obj.betas.b25;  % Multiply v_r with b25
            term1 = v_r .* term1;         % Element-wise multiplication

            term2 = v_r * obj.betas.b1;  
            term2 = term2 .* v_shift;  

            % --- Linear terms (growth/sinking) ---
            term3 = obj.linear * v_pos;  

            % --- Disaggregation terms (disabled for this version) ---
            term4 = 0 * v_pos;  

            % --- Combine all terms ---
            dvdt = (term1 + term2)' + term3 + term4;

            % --- Legacy disaggregation loop correction ---
            c3 = obj.config.c3;
            c4 = obj.config.c4;  
            for isec = 2:(n_sections-1)
                dvdt(isec) = dvdt(isec) - ...
                    c3 * c4^isec * (v_pos(isec) - c4 * v_pos(isec+1));
            end

            % Add NPP source term (Primary Production) ========
            if isprop(obj.config, 'use_NPP') && obj.config.use_NPP
                n = numel(v);
                s = zeros(n,1);
            
                % --- NEW block: read amplitude parameter if available ---
                ampNPP = 0.5; % default 50% amplitude
                if isprop(obj.config,'NPP_amp') && ~isempty(obj.config.NPP_amp)
                    ampNPP = obj.config.NPP_amp;
                end
                % -------------------------------------------------------

                % Time-varying productivity forcing 
                % Optionally modulate NPP as a function of time (sinusoidal seasonal driver)
                if isprop(obj.config,'use_NPP') && obj.config.use_NPP && ...
                        strcmpi(obj.config.NPP_profile,'sine')
                    % Sinusoidal variation around mean rate
                    rate = obj.config.NPP_rate * (1 + ampNPP*sin(2*pi*t/obj.config.t_final)); % updated
                else
                    rate = obj.config.NPP_rate; % fallback default
                end
                
                % Determine rate based on NPP profile type
                switch obj.config.NPP_profile
                    case 'constant'
                        rate = obj.config.NPP_rate;
                
                    case 'sine'
                        % Smooth daily-scale NPP oscillation over full experiment
                        T = obj.config.t_final;      % total days in experiment (e.g., 30)
                        amp = ampNPP;                % updated amplitude
                        rate = obj.config.NPP_rate * (1 + amp * sin(2*pi*t/T));
                
                    case 'step'
                        if t < obj.config.NPP_t_step
                            rate = obj.config.NPP_rate;
                        else
                            rate = obj.config.NPP_rate_after;
                        end
                
                    case 'pulse'
                        if abs(t - obj.config.NPP_t_step) < 0.5
                            rate = obj.config.NPP_rate;
                        else
                            rate = 0;
                        end
                
                    otherwise
                        rate = 0;
                end

                % --- Target section for source injection (default: 1) ---
                sec = 1;
                if isprop(obj.config, 'NPP_section')
                    sec = obj.config.NPP_section;
                end
                if sec >= 1 && sec <= n
                    s(sec) = rate;
                end

                % --- Add NPP source to ODE derivative ---
                dvdt = dvdt + s;
            end
        end

        function J = jacobian(obj, t, v)
            % JACOBIAN Evaluate analytical Jacobian
            % t = time (unused but required by ODE solver)
            % v = concentration vector
            % Returns: J = Jacobian matrix (derivatives of right-hand side with respect to variables)

            n_sections = length(v);
            v_r = v';

            % Build matrices for efficient computation
            v_mat = v_r(ones(1, n_sections), :);
            v_shift = [zeros(n_sections, 1), v_mat(:, 1:end-1)];

            % Term 1: diagonal (self-collisions)
            term1 = v_r * obj.betas.b25;
            term1 = diag(term1) + diag(v_r) .* obj.betas.b25;

            % Term 2: subdiagonal (collisions with smaller bin)
            term2a = v_r * obj.betas.b1;
            term2a = diag(term2a(2:end), -1);

            term2b = diag(obj.betas.b1, 1);
            term2b = term2b' .* v_r(1:end-1);
            term2b = diag(term2b, -1);

            term2c = diag(v_r(2:end), -1) .* obj.betas.b25';
            term2 = term2a + term2b + term2c;

            % Term 3: off-diagonal
            term3a = obj.betas.b1 .* v_shift;
            term3b = obj.betas.b25 .* v_mat;
            term3 = (term3a + term3b)';
            term3 = triu(term3, 2) + tril(term3, -1);

            % Assemble Jacobian
            J = term1 + term2 + term3 + obj.linear;

            % Include disaggregation matrices
            if n_sections > 2
                J = J - obj.disaggMinus + obj.disaggPlus;
            end
        end

        function [term1, term2, term3, term4, term5] = rateTerms(obj, v)
            % RATETERMS Evaluate individual rate terms for diagnostics
            n_sections = length(v);

            v_r = v';
            v_shift = [0, v_r(1:n_sections-1)];

            % Coagulation terms
            term1 = v_r * obj.betas.b25;
            term1 = v_r .* term1;

            term2 = v_r * obj.betas.b1;
            term2 = term2 .* v_shift;

            % Linear terms
            term3 = obj.linear * v;

            % Disaggregation terms
            term4 = -obj.disaggMinus * v;
            term5 = obj.disaggPlus * v;
        end

        function validate(obj)
            % VALIDATE Validate RHS configuration
            if isempty(obj.betas) || isempty(obj.linear)
                error('RHS not properly initialized');
            end

            n_sections = obj.config.n_sections;
            if size(obj.linear, 1) ~= n_sections
                error('Linear matrix dimension mismatch');
            end

            if obj.betas.getNumSections() ~= n_sections
                error('Beta matrices dimension mismatch');
            end
        end
    end
end