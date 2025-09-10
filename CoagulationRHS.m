classdef CoagulationRHS < handle
    %COAGULATIONRHS ODE right-hand side for coagulation equations

    properties
        betas;          % BetaMatrices object
        linear;         % Linear matrix (growth - sinking)
        disaggMinus;    % Disaggregation loss matrix
        disaggPlus;     % Disaggregation gain matrix
        config;         % SimulationConfig object
    end

    methods
        function obj = CoagulationRHS(betas, linear, disaggMinus, disaggPlus, config)
            %COAGULATIONRHS Constructor
            obj.betas = betas;
            obj.linear = linear;
            obj.disaggMinus = disaggMinus;
            obj.disaggPlus = disaggPlus;
            obj.config = config;
        end

        function dvdt = evaluate(obj, t, v)
            %EVALUATE Evaluate ODE right-hand side
            % t = time (unused but required by ODE solver)
            % v = concentration vector
            % Returns: dvdt = rate of change

            n_sections = length(v);

            % Match legacy behavior: clamp negative values to eps only for calculations
            v_pos = max(v, eps);

            % Convert to row vector for matrix operations - use clamped values
            v_r = v_pos';
            v_shift = [0, v_r(1:n_sections-1)];

            % Coagulation terms
            term1 = v_r * obj.betas.b25;
            term1 = v_r .* term1;

            term2 = v_r * obj.betas.b1;
            term2 = term2 .* v_shift;

            % Linear terms - use clamped values
            term3 = obj.linear * v_pos;

            % Disaggregation terms (disabled to match legacy RHS)
            term4 = 0*v_pos;

            % Combine all terms
            dvdt = (term1 + term2)' + term3 + term4;
            % Add legacy disaggregation loop contribution - use clamped values
            c3 = obj.config.c3;
            c4 = obj.config.c4;
            for isec = 2:(n_sections-1)
                dvdt(isec) = dvdt(isec) - c3 * c4^isec * (v_pos(isec) - c4 * v_pos(isec+1));
            end
        end

        function J = jacobian(obj, t, v)
            %JACOBIAN Evaluate analytical Jacobian
            % t = time (unused but required by ODE solver)
            % v = concentration vector
            % Returns: J = Jacobian matrix

            n_sections = length(v);
            v_r = v';

            % Build matrices for efficient computation
            v_mat = v_r(ones(1, n_sections), :);
            v_shift = [zeros(n_sections, 1), v_mat(:, 1:end-1)];

            % Term 1: df_i/dy_i terms
            term1 = v_r * obj.betas.b25;
            term1 = diag(term1) + diag(v_r) .* obj.betas.b25;

            % Term 2: df_i/dy_{i-1} terms
            term2a = v_r * obj.betas.b1;
            term2a = diag(term2a(2:end), -1);

            term2b = diag(obj.betas.b1, 1);
            term2b = term2b' .* v_r(1:end-1);
            term2b = diag(term2b, -1);

            term2c = diag(v_r(2:end), -1) .* obj.betas.b25';
            term2 = term2a + term2b + term2c;

            % Term 3: off-diagonal terms
            term3a = obj.betas.b1 .* v_shift;
            term3b = obj.betas.b25 .* v_mat;
            term3 = (term3a + term3b)';
            term3 = triu(term3, 2) + tril(term3, -1);

            % Assemble Jacobian
            % Match legacy Jacobian exactly: term1 + term2 + term3 + lin_term - disagg_minus + disagg_plus
            J = term1 + term2 + term3 + obj.linear;

            % Add disaggregation terms exactly as in legacy
            if n_sections > 2
                J = J - obj.disaggMinus + obj.disaggPlus;
            end
        end

        function [term1, term2, term3, term4, term5] = rateTerms(obj, v)
            %RATETERMS Evaluate individual rate terms for diagnostics
            % v = concentration vector
            % Returns: individual terms for analysis

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
            %VALIDATE Validate RHS configuration
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
