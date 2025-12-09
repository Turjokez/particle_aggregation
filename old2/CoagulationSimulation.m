classdef CoagulationSimulation < handle
    %COAGULATIONSIMULATION Main simulation controller for coagulation system

    properties
        config;         % SimulationConfig object
        grid;           % DerivedGrid object
        assembler;      % BetaAssembler object
        operators;      % Struct with linear operators
        rhs;            % CoagulationRHS object
        solver;         % ODESolver object
        result;         % Simulation results
    end

    methods
        function obj = CoagulationSimulation(varargin)
            %COAGULATIONSIMULATION Constructor
            % Can accept SimulationConfig object or parameter-value pairs

            if nargin == 0
                % Use default configuration
                obj.config = SimulationConfig();
            elseif nargin == 1 && isa(varargin{1}, 'SimulationConfig')
                % Use provided configuration
                obj.config = varargin{1};
            else
                % Create configuration from parameter-value pairs
                obj.config = SimulationConfig(varargin{:});
            end

            % Validate configuration
            obj.config.validate();

            % Initialize components
            obj.initializeComponents();
        end

        function initializeComponents(obj)
            %INITIALIZECOMPONENTS Initialize all simulation components

            % Create derived grid
            obj.grid = obj.config.derive();

            % Create beta assembler
            obj.assembler = BetaAssembler(obj.config, obj.grid);

            % Build linear operators
            obj.operators = struct();
            obj.operators.growth = LinearProcessBuilder.growthMatrix(obj.config, obj.grid);
            obj.operators.sink_loss = LinearProcessBuilder.sinkingMatrix(obj.config, obj.grid);
            [obj.operators.disagg_minus, obj.operators.disagg_plus] = ...
                LinearProcessBuilder.disaggregationMatrices(obj.config);
            obj.operators.linear = LinearProcessBuilder.linearMatrix(obj.config, obj.grid);

            % Create ODE solver
            obj.solver = ODESolver();

            % Initialize result structure
            obj.result = struct();
        end

        function result = run(obj, varargin)
            %RUN Execute the complete coagulation simulation
            % Optional: custom time span, initial conditions, solver options
            % Returns: result struct with simulation data

            fprintf('Starting coagulation simulation...\n');

            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'tspan', [], @isnumeric);
            addParameter(p, 'v0', [], @isnumeric);
            addParameter(p, 'solver_options', [], @isstruct);
            parse(p, varargin{:});

            % Set up time span
            if isempty(p.Results.tspan)
                tspan = obj.config.t_init:obj.config.delta_t:obj.config.t_final-1;
            else
                tspan = p.Results.tspan;
            end

            % Set up initial conditions
            if isempty(p.Results.v0)
                v0 = InitialSpectrumBuilder.initialSpectrum(obj.config, obj.grid);
            else
                v0 = p.Results.v0;
            end

            % Compute beta matrices
            fprintf('Computing coagulation kernels...\n');
            b_brown = obj.assembler.computeFor('KernelBrown');
            b_shear = obj.assembler.computeFor('KernelCurSh');
            b_ds = obj.assembler.computeFor('KernelCurDS');

            % Combine and scale beta matrices
            betas = obj.assembler.combineAndScale(b_brown, b_shear, b_ds);

            % Add betas to operators for diagnostics
            obj.operators.betas = betas;

            % Create RHS
            obj.rhs = CoagulationRHS(betas, obj.operators.linear, ...
                obj.operators.disagg_minus, obj.operators.disagg_plus, obj.config);

            % Validate RHS
            obj.rhs.validate();

            % Solve ODEs
            fprintf('Solving ODEs with disaggregation...\n');

            Y_current = v0(:);
            t = tspan(:);
            nsteps = length(t);
            Y = zeros(nsteps, length(Y_current));
            Y(1,:) = Y_current';
            
            for i = 2:nsteps
                % --- Coagulation step ---
                [~, Ytemp] = obj.solver.solve(obj.rhs, [t(i-1) t(i)], Y_current, p.Results.solver_options);
                Y_new = Ytemp(end,:)';  % last value from this step
            
                % =====================================================
                %  UPDATE: Time-varying epsilon(t) support
                % =====================================================
                eps_here = obj.config.epsilon; % default (constant)
                if isprop(obj.config,'epsilon_profile') && strcmpi(obj.config.epsilon_profile,'sine')
                    if isprop(obj.config,'epsilon_mean') && isprop(obj.config,'epsilon_amp') && ...
                       isprop(obj.config,'epsilon_period') && isprop(obj.config,'epsilon_phase')
                        eps_here = obj.config.epsilon_mean + obj.config.epsilon_amp * ...
                                   sin(2*pi*t(i)/obj.config.epsilon_period + obj.config.epsilon_phase);
                        eps_here = max(eps_here, 0); % prevent negative values
                    end
                end

                % === NEW ADDITION (Dynamic ε(t) scaling) ===
                % Rescale shear-driven coagulation dynamically
                if isfield(obj.operators,'betas')
                    try
                        beta_dynamic = obj.operators.betas; % copy base kernels
                        if isfield(beta_dynamic,'shear')
                            shear_scale = sqrt(eps_here / obj.config.epsilon_mean);
                            beta_dynamic.shear = beta_dynamic.shear * shear_scale;
                        end
                        if ismethod(obj.rhs,'updateKernel')
                            obj.rhs.updateKernel(beta_dynamic);
                        end
                    catch ME
                        warning('Dynamic kernel update skipped: %s', ME.message);
                    end
                end
                % =====================================================
            
                % --- Apply Disaggregation ---
                % === Nonlinear disaggregation feedback ===
                n_exp = 0.8; % exponent (~0.4–0.5 typical)
                Y_new = Disaggregation.apply(Y_new, obj.grid, eps_here * (eps_here / obj.config.epsilon_mean)^n_exp);
            
                % --- Store & advance ---
                Y(i,:) = Y_new';
                Y_current = Y_new;
            
                if mod(i,10)==0
                    fprintf('Step %d/%d → Disaggregation applied (ε=%.1e)\n', i, nsteps, eps_here);
                end
            end
            fprintf('Mass conserved check → Initial: %.3e | Final: %.3e\n', sum(Y(1,:)), sum(Y(end,:)));
            
            % Store results
            obj.result.time = t;
            obj.result.concentrations = Y;
            obj.result.initial_conditions = v0;
            obj.result.betas = betas;
            obj.result.operators = obj.operators;

            % Compute diagnostics
            fprintf('Computing diagnostics...\n');
            obj.result.diagnostics = obj.computeDiagnostics(t, Y);

            % Compute output data
            fprintf('Computing output data...\n');
            obj.result.output_data = OutputGenerator.spectraAndFluxes(t, Y, obj.grid, obj.config);

            fprintf('Simulation completed successfully.\n');
            result = obj.result;
        end

        function diagnostics = computeDiagnostics(obj, t, Y)
            %COMPUTEDIAGNOSTICS Compute simulation diagnostics
            % t = time vector
            % Y = concentration matrix
            % Returns: diagnostics struct

            diagnostics = struct();

            % Mass balance analysis
            [diagnostics.sectional_gains, diagnostics.sectional_losses] = ...
                MassBalanceAnalyzer.sectional(Y, obj.operators);

            [diagnostics.total_gains, diagnostics.total_losses] = ...
                MassBalanceAnalyzer.total(Y, obj.operators);

            % Rate terms over time
            diagnostics.rate_terms = obj.computeRateTermsOverTime(Y);

            % Mass conservation check
            diagnostics.mass_conservation = obj.checkMassConservation(Y);
        end

        function rate_terms = computeRateTermsOverTime(obj, Y)
            %COMPUTERATETERMSOVERTIME Compute rate terms for each time point
            % Y = concentration matrix
            % Returns: struct with rate terms over time

            n_times = size(Y, 1);
            n_sections = size(Y, 2);

            rate_terms = struct();
            rate_terms.term1 = zeros(n_times, n_sections);
            rate_terms.term2 = zeros(n_times, n_sections);
            rate_terms.term3 = zeros(n_times, n_sections);
            rate_terms.term4 = zeros(n_times, n_sections);
            rate_terms.term5 = zeros(n_times, n_sections);

            for i = 1:n_times
                [term1, term2, term3, term4, term5] = obj.rhs.rateTerms(Y(i, :)');
                rate_terms.term1(i, :) = term1';
                rate_terms.term2(i, :) = term2';
                rate_terms.term3(i, :) = term3';
                rate_terms.term4(i, :) = term4';
                rate_terms.term5(i, :) = term5';
            end
        end

        function conservation = checkMassConservation(obj, Y)
            %CHECKMASSCONSERVATION Check mass conservation
            % Y = concentration matrix
            % Returns: conservation struct with checks

            conservation = struct();

            % Total mass over time
            conservation.total_mass = sum(Y, 2);

            % Mass change rate (only if we have more than one time point)
            if size(Y, 1) > 1
                conservation.mass_change_rate = diff(conservation.total_mass);

                % Relative mass change (only if we have more than one time point)
                if length(conservation.total_mass) > 1
                    conservation.relative_change = conservation.mass_change_rate ./ conservation.total_mass(1:end-1);
                else
                    conservation.relative_change = [];
                end
            else
                conservation.mass_change_rate = [];
                conservation.relative_change = [];
            end

            % Conservation flag (mass should generally decrease due to settling)
            conservation.is_conserved = all(conservation.total_mass >= 0);

            % Settling loss estimate (per time: total across sections)
            if isfield(obj.operators, 'sink_loss')
                try
                    sink_diag = diag(obj.operators.sink_loss);
                    % Ensure column vector and matching section dimension
                    if isrow(sink_diag)
                        sink_diag = sink_diag';
                    end
                    if length(sink_diag) == size(Y, 2)
                        % Matrix product avoids element-wise broadcasting issues
                        conservation.settling_loss = Y * sink_diag;
                    else
                        conservation.settling_loss = [];
                    end
                catch ME
                    warning('Could not calculate settling loss: %s', ME.message);
                    conservation.settling_loss = [];
                end
            end
        end

        function generateOutputs(obj, plot_flag)
            %GENERATEOUTPUTS Generate all outputs and visualizations
            % plot_flag = whether to generate plots (default: true)

            if nargin < 2
                plot_flag = true;
            end

            if isempty(obj.result)
                error('No simulation results available. Run simulation first.');
            end

            % Generate plots if requested
            if plot_flag
                fprintf('Generating plots...\n');

                % Combine sectional gains and losses like legacy version
                % Legacy: sectional_gains = sec_gains.coag + sec_gains.growth;
                %         sectional_loss  = sec_losses.coag + sec_losses.settl + sec_losses.growth;
                combined_sectional_gains = obj.result.diagnostics.sectional_gains.coag + ...
                    obj.result.diagnostics.sectional_gains.growth;
                combined_sectional_losses = obj.result.diagnostics.sectional_losses.coag + ...
                    obj.result.diagnostics.sectional_losses.settl + ...
                    obj.result.diagnostics.sectional_losses.growth;

                OutputGenerator.plotAll(obj.result.time, obj.result.concentrations, ...
                    obj.result.output_data, obj.result.diagnostics.total_gains, ...
                    obj.result.diagnostics.total_losses, obj.config, ...
                    combined_sectional_gains, combined_sectional_losses, ...
                    obj.result.betas);
            end

            % Display diagnostics summary
            obj.displayDiagnosticsSummary();

            % Export data
            % obj.exportResults();
        end

        function displayDiagnosticsSummary(obj)
            %DISPLAYDIAGNOSTICSSUMMARY Display summary of simulation diagnostics

            fprintf('\n=== Simulation Diagnostics Summary ===\n');

            % Time and size info
            fprintf('Simulation time: %.2f to %.2f days\n', ...
                obj.result.time(1), obj.result.time(end));
            fprintf('Number of time points: %d\n', length(obj.result.time));
            fprintf('Number of sections: %d\n', size(obj.result.concentrations, 2));

            % Mass balance summary
            MassBalanceAnalyzer.displayBalanceSummary(...
                obj.result.diagnostics.sectional_gains, ...
                obj.result.diagnostics.sectional_losses, ...
                obj.result.time);

            % Conservation check
            if obj.result.diagnostics.mass_conservation.is_conserved
                fprintf('Mass conservation: PASSED\n');
            else
                fprintf('Mass conservation: FAILED\n');
            end

            % Beta matrices summary
            if isfield(obj.result, 'betas')
                obj.result.betas.displaySummary();
            end
        end

        function exportResults(obj, filename)
            %EXPORTRESULTS Export simulation results to file
            % filename = optional output filename

            if nargin < 2
                filename = sprintf('coagulation_simulation_%s.mat', ...
                    datestr(now, 'yyyymmdd_HHMMSS'));
            end

            % Export output data
            OutputGenerator.exportData(obj.result.output_data, filename);

            % Export full results
            results_filename = strrep(filename, '.mat', '_full.mat');
            save(results_filename, 'obj');
            fprintf('Full simulation results saved to: %s\n', results_filename);
        end

        function enableTracer(obj)
            %ENABLETRACER Enable tracer integration (future implementation)
            warning('Tracer integration not yet implemented');
        end
    end
end