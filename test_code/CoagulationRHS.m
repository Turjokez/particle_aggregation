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
            
            % 1. Physics Scaling (Epsilon & Alpha)
            [coag_scale, ~] = obj.getScalingFactors(t);
            
            % 2. Prepare State
            v    = v(:);                  % Column vector
            v_pos = max(v, eps);          % Positivity guard
            
            Ns = obj.config.n_sections;
            
            if obj.config.use_column
                % === 1-D COLUMN MODE ===
                Nz = floor(obj.config.z_max / obj.config.dz);
                V_grid = reshape(v_pos, Ns, Nz).'; % [Nz x Ns]
                dCoag_grid = zeros(Nz, Ns);
                b25 = obj.betas.b25;
                b1  = obj.betas.b1;
                
                for k = 1:Nz
                    row_v       = V_grid(k, :);          % 1 x Ns
                    row_v_shift = [0, row_v(1:end-1)];   % shifted spectrum
                    term1 = row_v .* (row_v * b25);
                    term2 = (row_v * b1) .* row_v_shift;
                    dCoag_grid(k, :) = term1 + term2;
                end
                
                dCoag_grid = dCoag_grid * coag_scale;
                dvdt_coag  = dCoag_grid.'; 
                dvdt_coag  = dvdt_coag(:);
                
                % Add Linear Transport (Sinking + Growth)
                dvdt = dvdt_coag + (obj.linear * v_pos);
                
            else
                % === 0-D SLAB MODE ===
                v_r     = v_pos.'; 
                v_shift = [0, v_r(1:end-1)];
                term1   = v_r .* (v_r * obj.betas.b25);
                term2   = (v_r * obj.betas.b1) .* v_shift;
                dvdt_coag = (term1 + term2).' * coag_scale;
                
                dvdt = dvdt_coag + (obj.linear * v_pos);
                
                if ~obj.config.disagg_use_nonlinear
                    dvdt = dvdt - (obj.disaggMinus * v_pos) + (obj.disaggPlus * v_pos);
                end
            end
            
            % 3. NPP Source
            if obj.config.use_NPP
                source_vec = obj.getNPPSource(t, Ns);
                dvdt = dvdt + source_vec;
            end
            
            % === NEW: ATTENUATION / GRAZING TERM =========================
            % Old behaviour: uniform mu * v_pos
            % New: same for 0-D; depth-dependent mu(z) for 1-D column.
            if isprop(obj.config, 'attenuation_rate') && obj.config.attenuation_rate > 0
                mu0 = obj.config.attenuation_rate;   % base rate [d^-1]
                
                if obj.config.use_column
                    % Depth dependence: mu(z) = mu0 * (1 + lambda * depth)
                    Nz = floor(obj.config.z_max / obj.config.dz);
                    
                    if isprop(obj.config, 'attenuation_depth_factor')
                        lambda = obj.config.attenuation_depth_factor;
                    else
                        lambda = 0.0;
                    end
                    
                    depth_vec = ((0:Nz-1)' + 0.5) * obj.config.dz;  % mid-layer depth [m]
                    mu_z      = mu0 * (1 + lambda * depth_vec);      % [Nz x 1]
                    
                    % Expand to all size bins
                    mu_big = kron(mu_z, ones(Ns,1));                % [Ns*Nz x 1]
                    dvdt   = dvdt - mu_big .* v_pos;
                else
                    % 0-D slab: uniform attenuation
                    dvdt = dvdt - mu0 * v_pos;
                end
            end
            
            % 4. Guard
            if any(~isfinite(dvdt))
                dvdt(~isfinite(dvdt)) = 0;
            end
        end
        
        function J = jacobian(obj, t, v) %#ok<INUSD>
            J = obj.linear;
        end
        
        function J_loc = computeLocalJacobian(obj, v, scale), J_loc=[]; end %#ok<INUSD>
        
        function [coag_scale, eps_here] = getScalingFactors(obj, t)
            if strcmpi(obj.config.epsilon_profile,'observed') && ~isempty(obj.config.epsilon_series)
                if ~isempty(obj.config.epsilon_time)
                    eps_here = interp1(obj.config.epsilon_time(:), obj.config.epsilon_series(:), ...
                                       t, 'linear', 'extrap');
                else
                    eps_here = obj.config.epsilon;
                end
            else
                eps_here = obj.config.epsilon; 
            end
            eps_here = max(eps_here, 1e-14);
            ratio    = eps_here / obj.config.epsilon_ref;
            alpha_eff = obj.config.alpha_base * (ratio ^ -obj.config.p_alpha);
            alpha_eff = min(max(alpha_eff, 1e-3), 50.0);
            coag_scale = alpha_eff; 
        end
        
        function s = getNPPSource(obj, t, Ns) %#ok<INUSD>
            rate = obj.config.NPP_rate; 
            if obj.config.use_column
                Nz = floor(obj.config.z_max / obj.config.dz);
                s_grid = zeros(Nz, Ns);
                s_grid(1, 1) = rate; 
                s = s_grid.'; s = s(:);
            else
                s = zeros(Ns, 1); 
                s(1) = rate;
            end
        end
        
        function validate(obj), end
        function [term1, term2, term3, term4, term5] = rateTerms(obj, v, t)
            n = length(v); term1=zeros(n,1); term2=term1; term3=term1; term4=term1; term5=term1;
        end
    end
end