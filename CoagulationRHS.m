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
            obj.betas = betas;  % Set beta matrices
            obj.linear = linear;  % Set linear matrix
            obj.disaggMinus = disaggMinus;  % Set disaggregation loss matrix
            obj.disaggPlus = disaggPlus;  % Set disaggregation gain matrix
            obj.config = config;  % Set configuration object
        end

        function dvdt = evaluate(obj, t, v)
            % EVALUATE Evaluate ODE right-hand side
            % t = time (unused but required by ODE solver)
            % v = concentration vector
            % Returns: dvdt = rate of change (the right-hand side of ODE)

            n_sections = length(v);  % Number of sections (length of concentration vector)

            % Match legacy behavior: clamp negative values to eps only for calculations
            v_pos = max(v, eps);  % Ensure that concentrations are non-negative by replacing negatives with a small value

            v_r = v_pos';  % Convert the concentration vector to a row vector for matrix operations
            v_shift = [0, v_r(1:n_sections-1)];  % Shift the concentration vector to align with previous section

            % Coagulation terms
            term1 = v_r * obj.betas.b25;  % Multiply v_r with b25 (coagulation kernel matrix)
            term1 = v_r .* term1;  % Element-wise multiplication

            term2 = v_r * obj.betas.b1;  
            term2 = term2 .* v_shift;  

            % Linear terms - growth and sinking
            term3 = obj.linear * v_pos;  

            % Disaggregation terms (disabled to match legacy RHS)
            term4 = 0 * v_pos;  % Set term4 to zero as disaggregation terms are disabled in this version

            % Combine all terms to calculate the rate of change
            dvdt = (term1 + term2)' + term3 + term4;

            % Add legacy disaggregation loop contribution - adjust for negative and positive shifts
            c3 = obj.config.c3;  % Coefficient c3 from config
            c4 = obj.config.c4;  
            for isec = 2:(n_sections-1)  % Loop through all sections except the first and last
                dvdt(isec) = dvdt(isec) - c3 * c4^isec * (v_pos(isec) - c4 * v_pos(isec+1));  % Apply disaggregation correction
            end
        end

        function J = jacobian(obj, t, v)
            % JACOBIAN Evaluate analytical Jacobian
            % t = time (unused but required by ODE solver)
            % v = concentration vector
            % Returns: J = Jacobian matrix (derivatives of right-hand side with respect to variables)

            n_sections = length(v);  % Number of sections (length of concentration vector)
            v_r = v';  % Convert concentration vector to a row vector

            % Build matrices for efficient computation
            v_mat = v_r(ones(1, n_sections), :);  % Replicate the row vector to create a matrix
            v_shift = [zeros(n_sections, 1), v_mat(:, 1:end-1)];  % Shift the concentration vector

            % Term 1: df_i/dy_i terms (diagonal elements of Jacobian)
            term1 = v_r * obj.betas.b25;  % Multiply v_r with b25 (coagulation kernel matrix)
            term1 = diag(term1) + diag(v_r) .* obj.betas.b25;  % Diagonal terms: self-collisions

            % Term 2: df_i/dy_{i-1} terms (subdiagonal elements of Jacobian)
            term2a = v_r * obj.betas.b1;  % Multiply v_r with b1 (coagulation kernel matrix)
            term2a = diag(term2a(2:end), -1);  % Shift and create diagonal matrix

            term2b = diag(obj.betas.b1, 1);  % Create diagonal matrix from b1 for next section
            term2b = term2b' .* v_r(1:end-1);  % Element-wise multiplication with shifted concentration vector
            term2b = diag(term2b, -1);  % Create diagonal matrix from shifted b1 terms

            term2c = diag(v_r(2:end), -1) .* obj.betas.b25';  % Apply b25 matrix for inter-section collisions
            term2 = term2a + term2b + term2c;  % Combine all terms for the subdiagonal elements

            % Term 3: Off-diagonal terms (terms where j ≠ i or j ≠ i-1)
            term3a = obj.betas.b1 .* v_shift;  % Element-wise multiplication for off-diagonal terms
            term3b = obj.betas.b25 .* v_mat;  % Element-wise multiplication for off-diagonal terms
            term3 = (term3a + term3b)';  % Combine and transpose for proper dimensions
            term3 = triu(term3, 2) + tril(term3, -1);  % Create upper and lower triangular matrices for off-diagonal terms

            % Assemble the Jacobian matrix
            J = term1 + term2 + term3 + obj.linear;  % Sum the individual terms to get the Jacobian matrix

            % Add disaggregation terms (if available) exactly as in legacy
            if n_sections > 2  % If more than two sections, include disaggregation terms
                J = J - obj.disaggMinus + obj.disaggPlus;  % Add disaggregation loss and gain matrices
            end
        end

        function [term1, term2, term3, term4, term5] = rateTerms(obj, v)
            % RATETERMS Evaluate individual rate terms for diagnostics
            % v = concentration vector
            % Returns: individual terms for analysis

            n_sections = length(v);  % Number of sections (length of concentration vector)

            v_r = v';  % Convert concentration vector to row vector
            v_shift = [0, v_r(1:n_sections-1)];  % Shift the concentration vector

            % Coagulation terms
            term1 = v_r * obj.betas.b25;  % Multiply v_r with b25 (coagulation kernel matrix)
            term1 = v_r .* term1;  % Element-wise multiplication

            term2 = v_r * obj.betas.b1;  % Multiply v_r with b1 (coagulation kernel matrix)
            term2 = term2 .* v_shift;  % Element-wise multiplication with shifted concentration vector

            % Linear terms - growth and sinking
            term3 = obj.linear * v;  % Linear term applied to concentration vector

            % Disaggregation terms - loss and gain
            term4 = -obj.disaggMinus * v;  % Loss due to disaggregation
            term5 = obj.disaggPlus * v;  % Gain due to disaggregation
        end

        function validate(obj)
            % VALIDATE Validate RHS configuration
            % Ensures that all matrices and configurations are properly initialized

            if isempty(obj.betas) || isempty(obj.linear)  % Check if required matrices are initialized
                error('RHS not properly initialized');
            end

            n_sections = obj.config.n_sections;  % Get the number of sections from the config
            if size(obj.linear, 1) ~= n_sections  % Check if the linear matrix has the correct dimensions
                error('Linear matrix dimension mismatch');
            end

            if obj.betas.getNumSections() ~= n_sections  % Check if beta matrices have the correct dimensions
                error('Beta matrices dimension mismatch');
            end
        end
    end
end
