classdef CoagulationRHS < handle
    % COAGULATIONRHS ODE right-hand side (Supports 0-D Slab and 1-D Column)
    
    properties
        betas;          % BetaMatrices object
        linear;         % Linear matrix (Growth + Sinking)
        disaggMinus;    % Linear Disagg (Legacy)
        disaggPlus;     % Linear Disagg (Legacy)
        config;         % SimulationConfig object
    end

    methods
        function obj = CoagulationRHS(betas, linear, disaggMinus, disaggPlus, config)
            obj.betas       = betas;
            obj.linear      = linear;
            obj.disaggMinus = disaggMinus;
            obj.disaggPlus  = disaggPlus;
            obj.config      = config;
        end

        function updateKernel(obj, newBetas)
            obj.betas = newBetas;
        end

        function dvdt = evaluate(obj, t, v)
            % EVALUATE Compute dy/dt
            % Handles both Slab (Ns x 1) and Column (Nz*Ns x 1)
            
            % 1. Physics Scaling (Epsilon & Alpha)
            [coag_scale, ~] = obj.getScalingFactors(t);
            
            % 2. Prepare State
            v = v(:); % Column vector
            v_pos = max(v, eps); % Positivity guard
            
            Ns = obj.config.n_sections;
            
            if obj.config.use_column
                % === 1-D COLUMN MODE ===
                Nz = floor(obj.config.z_max / obj.config.dz);
                
                % Reshape to [Layers x Sections]
                % We assume state vector is ordered: Layer1, Layer2...
                % Note: reshape fills columns first. If state is [L1_s1, L1_s2...; L2_s1...]
                % Then we reshape to [Ns, Nz] and transpose.
                V_grid = reshape(v_pos, Ns, Nz).'; % Result: [Nz x Ns]
                
                % Calculate Coagulation for each layer
                % (We loop because Beta matrices are constant, but Concentrations vary)
                dCoag_grid = zeros(Nz, Ns);
                
                b25 = obj.betas.b25;
                b1  = obj.betas.b1;
                
                % Pre-calculate shift indices for speed
                for k = 1:Nz
                    row_v = V_grid(k, :);     % 1 x Ns
                    
                    % Shifted vector (for Gain term)
                    row_v_shift = [0, row_v(1:end-1)];
                    
                    % Term 1: Gain/Loss from collisions (Matrix Vector)
                    % b25 is [Ns x Ns]. row_v is [1 x Ns].
                    % Formula: v .* (v * b25)
                    term1 = row_v .* (row_v * b25);
                    
                    % Term 2: Gain from smaller particles
                    % Formula: (v * b1) .* v_shift
                    term2 = (row_v * b1) .* row_v_shift;
                    
                    dCoag_grid(k, :) = term1 + term2;
                end
                
                % Apply Physics Scaling (turbulence effect on stickiness)
                dCoag_grid = dCoag_grid * coag_scale;
                
                % Flatten back to column vector
                dvdt_coag = dCoag_grid.'; % [Ns x Nz]
                dvdt_coag = dvdt_coag(:); % [Total x 1]
                
                % Add Linear Transport (Sinking + Growth)
                % obj.linear is the big sparse matrix from LinearProcessBuilder
                dvdt = dvdt_coag + (obj.linear * v_pos);
                
            else
                % === 0-D SLAB MODE (Legacy Optimized) ===
                v_r = v_pos.'; 
                v_shift = [0, v_r(1:end-1)];
                
                term1 = v_r .* (v_r * obj.betas.b25);
                term2 = (v_r * obj.betas.b1) .* v_shift;
                
                dvdt_coag = (term1 + term2).' * coag_scale;
                
                dvdt = dvdt_coag + (obj.linear * v_pos);
                
                % Legacy Linear Disagg (if not using nonlinear impulse)
                if ~obj.config.disagg_use_nonlinear
                     dvdt = dvdt - (obj.disaggMinus * v_pos) + (obj.disaggPlus * v_pos);
                end
            end
            
            % 3. NPP Source (Primary Production)
            % For Column: Inject into top layer(s) only
            if obj.config.use_NPP
                 source_vec = obj.getNPPSource(t, Ns);
                 dvdt = dvdt + source_vec;
            end
            
            % 4. NaN Guard
            if any(~isfinite(dvdt))
                dvdt(~isfinite(dvdt)) = 0;
            end
        end

        % =================================================================
        % === ANALYTICAL JACOBIAN (Required for ODE15s stability) =========
        % =================================================================
        function J = jacobian(obj, t, v)
            % Returns the Jacobian d(f)/d(v)
            [coag_scale, ~] = obj.getScalingFactors(t);
            v = v(:);
            v_pos = max(v, eps);
            Ns = obj.config.n_sections;
            
            if obj.config.use_column
                % === 1-D Jacobian (Block Diagonal Coag + Global Linear) ===
                Nz = floor(obj.config.z_max / obj.config.dz);
                
                % We need to build a block-diagonal matrix where each block is the 
                % 0-D coagulation Jacobian for that layer.
                
                % 1. Reshape state
                V_grid = reshape(v_pos, Ns, Nz).'; % [Nz x Ns]
                
                % 2. Pre-allocate sparse triplet arrays
                max_entries = Nz * Ns * Ns;
                rows = zeros(max_entries, 1);
                cols = zeros(max_entries, 1);
                vals = zeros(max_entries, 1);
                
                idx_ptr = 0;
                
                % 3. Loop over layers to build local Jacobians
                for k = 1:Nz
                    v_local = V_grid(k, :).'; % [Ns x 1]
                    
                    % Calculate Local Jacobian J_local [Ns x Ns]
                    J_local = obj.computeLocalJacobian(v_local, coag_scale);
                    
                    % Store in triplet lists
                    [r, c, val] = find(J_local);
                    n_vals = length(val);
                    
                    offset = (k-1) * Ns;
                    
                    rows(idx_ptr+1 : idx_ptr+n_vals) = r + offset;
                    cols(idx_ptr+1 : idx_ptr+n_vals) = c + offset;
                    vals(idx_ptr+1 : idx_ptr+n_vals) = val;
                    
                    idx_ptr = idx_ptr + n_vals;
                end
                
                % 4. Create Sparse Block Matrix
                J_coag = sparse(rows(1:idx_ptr), cols(1:idx_ptr), vals(1:idx_ptr), ...
                                Nz*Ns, Nz*Ns);
                            
                % 5. Add Linear Transport (Sinking)
                J = J_coag + obj.linear;
                
            else
                % === 0-D Jacobian ===
                J_coag = obj.computeLocalJacobian(v_pos, coag_scale);
                J = J_coag + obj.linear;
                
                % Legacy Disagg
                if ~obj.config.disagg_use_nonlinear
                     J = J - obj.disaggMinus + obj.disaggPlus;
                end
            end
        end

        % Helper: Compute 0-D Coagulation Jacobian for one layer
        function J_loc = computeLocalJacobian(obj, v, scale)
            % v is column [Ns x 1]
            Ns = length(v);
            v_r = v.'; 
            
            % Standard Coagulation Jacobian Construction
            
            v_mat   = v_r(ones(1, Ns), :);         % Ns x Ns row-repeated
            v_shift = [zeros(Ns, 1), v_mat(:, 1:end-1)];
            
            % Term 1
            t1_vec = v_r * obj.betas.b25;
            term1 = diag(t1_vec) + diag(v_r) .* obj.betas.b25;
            
            % Term 2
            t2a = v_r * obj.betas.b1; 
            term2a = diag(t2a(2:end), -1);
            
            % === FIX: Corrected variable naming here ===
            term2b = diag(obj.betas.b1, 1);
            term2b = term2b.' .* v_r(1:end-1);
            term2b = diag(term2b, -1);
            
            t2c = diag(v_r(2:end), -1) .* obj.betas.b25.';
            
            term2 = term2a + term2b + t2c;
            
            % Term 3
            term3a = obj.betas.b1  .* v_shift;
            term3b = obj.betas.b25 .* v_mat;
            term3  = (term3a + term3b).';
            term3  = triu(term3, 2) + tril(term3, -1);
            
            J_loc = scale * (term1 + term2 + term3);
        end

        % Helper: Calculate Alpha and Shear Scaling from Epsilon
        function [coag_scale, eps_here] = getScalingFactors(obj, t)
            % 1. Get Epsilon
            if strcmpi(obj.config.epsilon_profile,'observed') && ~isempty(obj.config.epsilon_series)
                 if ~isempty(obj.config.epsilon_time)
                    eps_here = interp1(obj.config.epsilon_time(:), obj.config.epsilon_series(:), ...
                                       t, 'linear', 'extrap');
                 else
                     eps_here = obj.config.epsilon;
                 end
            else
                % Constant or Sine
                eps_here = obj.config.epsilon; 
                if strcmpi(obj.config.epsilon_profile,'sine')
                    eps_here = obj.config.epsilon_mean + obj.config.epsilon_amp * ...
                           sin(2*pi*t/obj.config.epsilon_period + obj.config.epsilon_phase);
                end
            end
            
            eps_here = max(eps_here, 1e-14);
            
            % 2. Stickiness Scaling
            ratio = eps_here / obj.config.epsilon_ref;
            alpha_eff = obj.config.alpha_base * (ratio ^ -obj.config.p_alpha);
            
            % Clip for stability
            alpha_eff = min(max(alpha_eff, 1e-3), 50.0);
            
            coag_scale = alpha_eff; 
        end
        
        % Helper: NPP Source Injection
        function s = getNPPSource(obj, t, Ns)
             rate = obj.config.NPP_rate; % Base rate
             
             if strcmpi(obj.config.NPP_profile, 'sine')
                 rate = rate * (1 + 0.5*sin(2*pi*t/30)); 
             end
             
             if obj.config.use_column
                 Nz = floor(obj.config.z_max / obj.config.dz);
                 s_grid = zeros(Nz, Ns);
                 s_grid(1, 1) = rate; % Top layer only
                 s = s_grid.'; 
                 s = s(:);
             else
                 s = zeros(Ns, 1);
                 s(1) = rate;
             end
        end

        function validate(obj)
            if isempty(obj.betas) || isempty(obj.linear)
                error('RHS not properly initialized');
            end
        end
        
        % Placeholder for rate terms
        function [term1, term2, term3, term4, term5] = rateTerms(obj, v, t)
             n = length(v);
             term1 = zeros(n,1); term2 = zeros(n,1); 
             term3 = zeros(n,1); term4 = zeros(n,1); term5 = zeros(n,1);
        end
    end
end