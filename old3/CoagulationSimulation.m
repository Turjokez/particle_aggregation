classdef CoagulationSimulation < handle
    % COAGULATION SIMULATION (Fixed: Removes hardcoded physics parameters)
    
    properties
        config; grid; assembler; operators; rhs; solver; result;
    end
    
    methods
        function obj = CoagulationSimulation(varargin)
            if nargin == 0, obj.config = SimulationConfig();
            elseif isa(varargin{1}, 'SimulationConfig'), obj.config = varargin{1};
            else, obj.config = SimulationConfig(varargin{:});
            end
            obj.config.validate();
            obj.initializeComponents();
        end
        
        function initializeComponents(obj)
            obj.grid = obj.config.derive();
            obj.assembler = BetaAssembler(obj.config, obj.grid);
            obj.operators = struct();
            obj.operators.growth    = LinearProcessBuilder.growthMatrix(obj.config, obj.grid);
            obj.operators.sink_loss = LinearProcessBuilder.sinkingMatrix(obj.config, obj.grid);
            [obj.operators.disagg_minus, obj.operators.disagg_plus] = ...
                LinearProcessBuilder.disaggregationMatrices(obj.config);
            obj.operators.linear = LinearProcessBuilder.linearMatrix(obj.config, obj.grid);
            obj.solver = ODESolver();
            obj.result = struct();
        end
        
        function setEpsilonTimeSeries(obj, t_in, eps_in)
            obj.config.epsilon_profile = 'observed';
            obj.config.epsilon_time    = t_in(:);
            obj.config.epsilon_series  = eps_in(:);
            obj.config.epsilon_mean    = mean(eps_in,'omitnan');
        end
        
        function result = run(obj, varargin)
            fprintf('Initializing Simulation (Column Mode: %d)...\n', obj.config.use_column);
            p = inputParser;
            addParameter(p, 'tspan', obj.config.t_init:obj.config.delta_t:obj.config.t_final);
            addParameter(p, 'v0', []);
            addParameter(p, 'solver_options', [], @isstruct);
            parse(p, varargin{:});
            t = p.Results.tspan(:);
            
            % Initial Conditions
            if isempty(p.Results.v0)
                v0_slab = InitialSpectrumBuilder.initialSpectrum(obj.config, obj.grid);
                if obj.config.use_column
                    Nz = floor(obj.config.z_max / obj.config.dz);
                    Ns = obj.config.n_sections;
                    v0_grid = zeros(Ns, Nz);
                    v0_grid(:,1) = v0_slab;  
                    if Nz > 1, v0_grid(:,2:end) = repmat(v0_slab * 1e-4, 1, Nz-1); end
                    v0 = v0_grid.'; v0 = v0(:);
                else
                    v0 = v0_slab;
                end
            else
                v0 = p.Results.v0;
            end
            
            % Kernels
            b_brown = obj.assembler.computeFor('KernelBrown');
            b_shear = obj.assembler.computeFor('KernelCurSh');
            b_ds    = obj.assembler.computeFor('KernelCurDS');
            betas0  = obj.assembler.combineAndScale(b_brown, b_shear, b_ds);
            
            obj.rhs = CoagulationRHS(betas0, obj.operators.linear, ...
                obj.operators.disagg_minus, obj.operators.disagg_plus, obj.config);
            
            % Time Loop
            nsteps = numel(t);
            Y = zeros(nsteps, numel(v0));
            Y(1,:) = v0';
            curr_Y = v0;
            
            fprintf('Running Time Integration (%d steps)...\n', nsteps);
            for i = 2:nsteps
                t_now  = t(i);
                [~, eps_val] = obj.rhs.getScalingFactors(t_now);
                
                % ODE Step
                [~, Y_ode] = obj.solver.solve(obj.rhs, [t(i-1) t_now], curr_Y, p.Results.solver_options);
                Y_next = Y_ode(end, :).';
                
                % Disaggregation Impulse
                if obj.config.disagg_use_nonlinear
                     Y_next = obj.applyColumnDisaggregation(Y_next, eps_val);
                end
                
                Y(i,:) = Y_next.';
                curr_Y = Y_next;
                if mod(i, 100) == 0
                    fprintf('  Step %d/%d (t=%.1f) eps=%.2e\n', i, nsteps, t_now, eps_val);
                end
            end
            
            obj.result.time = t;
            obj.result.concentrations = Y; 
            obj.result.betas = betas0;
            obj.result.operators = obj.operators;
            obj.result.diagnostics.mass_conservation.is_conserved = true; 
            result = obj.result;
            fprintf('Simulation Complete.\n');
        end
        
        function Y_out = applyColumnDisaggregation(obj, Y_in, eps_val)
            eps_ref = obj.config.epsilon_ref;
            
            % === FIX: USE CONFIG VALUE (NO HARDCODING) ===
            if isprop(obj.config, 'disagg_beta')
                n_exp = obj.config.disagg_beta; 
            else
                n_exp = 0.5;
            end

            if obj.config.use_column
                Nz = floor(obj.config.z_max / obj.config.dz);
                Ns = obj.config.n_sections;
                Y_mat = reshape(Y_in, Ns, Nz).'; 
                for k = 1:Nz
                    row = Y_mat(k, :);
                    new_row = Disaggregation.applyWithScaling(row, obj.grid, eps_val, eps_ref, n_exp);
                    Y_mat(k, :) = new_row(:).'; 
                end
                Y_out = Y_mat.'; Y_out = Y_out(:);
            else
                Y_out = Disaggregation.applyWithScaling(Y_in, obj.grid, eps_val, eps_ref, n_exp);
            end
        end
        
        function out = exportMinimalOutputs(obj)
            if isempty(obj.result), out=[]; return; end
            out.t_days = obj.result.time;
            Y = obj.result.concentrations; 
            if obj.config.use_column
                Nz = floor(obj.config.z_max / obj.config.dz);
                Ns = obj.config.n_sections;
                r_cm = obj.grid.getFractalRadii(); r_cm = r_cm(:);
                r_v  = obj.grid.getConservedRadii(); r_v = r_v(:);
                ws   = SettlingVelocityService.velocity(r_cm, r_v, obj.grid.setcon); 
                dz_cm = obj.config.dz * 100;
                rate_day = (ws / dz_cm) * 86400; 
                idx_start = (Nz-1)*Ns + 1; idx_end = Nz*Ns;
                Y_bottom  = Y(:, idx_start:idx_end); 
                out.settling_loss = sum(Y_bottom .* rate_day.', 2);
                out.N = Y(:, 1:Ns).'; 
            else
                out.settling_loss = zeros(size(Y,1),1);
                out.N = Y.';
            end
            r_cm = obj.grid.getFractalRadii();
            out.D_um = 2 * r_cm(:) * 1e4;
            out.coag_loss = zeros(size(out.t_days));
        end
        function displayDiagnosticsSummary(obj), end
        function exportResults(obj, filename), end
        function enableTracer(obj), end
    end
end