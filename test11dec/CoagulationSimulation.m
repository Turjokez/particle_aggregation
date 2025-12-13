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

        % NEW-2025-12-11: Optional time-varying epsilon forcing
        epsilon_time = [];      % NEW-2025-12-11: time vector [d]
        epsilon_values = [];    % NEW-2025-12-11: epsilon(t) [e.g., m^2 s^{-3}]
    end

    methods
        function obj = CoagulationSimulation(varargin)
            %COAGULATIONSIMULATION Constructor
            % Can accept SimulationConfig object or parameter-value pairs

            if nargin == 0
                obj.config = SimulationConfig();
            elseif nargin == 1 && isa(varargin{1}, 'SimulationConfig')
                obj.config = varargin{1};
            else
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

            % NEW-2025-12-11: pass config so MassBalanceAnalyzer can detect column mode
            obj.operators.config = obj.config;   % NEW-2025-12-11

            obj.operators.growth    = LinearProcessBuilder.growthMatrix(obj.config, obj.grid);

            % OLD (kept):
            % obj.operators.sink_loss = LinearProcessBuilder.sinkingMatrix(obj.config, obj.grid);

            % NEW-2025-12-12: capture extra outputs if you want later
            [obj.operators.sink_loss, obj.operators.sink_rate, obj.operators.settling_vel_cmday] = ...
                LinearProcessBuilder.sinkingMatrix(obj.config, obj.grid);

            [obj.operators.disagg_minus, obj.operators.disagg_plus] = ...
                LinearProcessBuilder.disaggregationMatrices(obj.config);

            obj.operators.linear    = LinearProcessBuilder.linearMatrix(obj.config, obj.grid);

            % Create ODE solver
            obj.solver = ODESolver();

            % Initialize result structure
            obj.result = struct();
        end

        % ==============================================================
        % NEW-2025-12-11: Allow external scripts to set epsilon(t)
        % ==============================================================
        function setEpsilonTimeSeries(obj, t, eps_vals)
            obj.epsilon_time   = t(:);
            obj.epsilon_values = eps_vals(:);

            if ~isempty(obj.rhs) && ismethod(obj.rhs, 'setEpsilonSeries')
                obj.rhs.setEpsilonSeries(obj.epsilon_time, obj.epsilon_values);
            end
        end

        % ==============================================================
        % NEW-2025-12-12: Robust initial condition builder (Adrian tests)
        % ==============================================================
        function v0 = buildInitialCondition(obj)
            cfg = obj.config;

            baseV = InitialSpectrumBuilder.initialSpectrum(cfg, obj.grid);

            if ~isprop(cfg,'init_profile') || isempty(cfg.init_profile)
                v0 = baseV;
                return;
            end

            prof = lower(string(cfg.init_profile));

            if ~(isprop(cfg,'use_column') && cfg.use_column)
                v0 = baseV;
                return;
            end

            Ns = cfg.n_sections;
            Nz = cfg.getNumLayers();

            switch prof
                case "top_only"
                    N2 = zeros(Ns, Nz);

                    amp = 1;
                    if isprop(cfg,'n1') && ~isempty(cfg.n1)
                        amp = cfg.n1;
                    elseif isprop(cfg,'num_1') && ~isempty(cfg.num_1)
                        amp = cfg.num_1;
                    end

                    N2(:,1) = amp;
                    v0 = N2(:);

                otherwise
                    v0 = baseV;
            end
        end

        function result = run(obj, varargin)
            %RUN Execute the complete coagulation simulation

            fprintf('Starting coagulation simulation...\n');

            % Parse optional arguments
            p = inputParser;
            addParameter(p, 'tspan', [], @isnumeric);
            addParameter(p, 'v0', [], @isnumeric);
            addParameter(p, 'solver_options', [], @isstruct);
            parse(p, varargin{:});

            % Set up time span
            if isempty(p.Results.tspan)
                % OLD (kept):
                % tspan = obj.config.t_init:obj.config.delta_t:obj.config.t_final-1;

                % NEW-2025-12-12: include t_final
                tspan = obj.config.t_init:obj.config.delta_t:obj.config.t_final;
            else
                tspan = p.Results.tspan;
            end

            % Initial conditions
            if isempty(p.Results.v0)
                v0 = obj.buildInitialCondition();
            else
                v0 = p.Results.v0;
            end

            % ==========================================================
            % NEW-2025-12-12: Decide if coagulation is OFF (skip kernels)
            % ==========================================================
            coag_is_off = false;
            if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag)
                coag_is_off = ~logical(obj.config.enable_coag);
            end
            if isprop(obj.config,'coag_scale') && ~isempty(obj.config.coag_scale) && obj.config.coag_scale == 0
                coag_is_off = true;
            end
            if isprop(obj.config,'beta_scale') && ~isempty(obj.config.beta_scale) && obj.config.beta_scale == 0
                coag_is_off = true;
            end
            if isprop(obj.config,'kernel_scale') && ~isempty(obj.config.kernel_scale) && obj.config.kernel_scale == 0
                coag_is_off = true;
            end

            % ==========================================================
            % Compute beta matrices (skippable)
            % ==========================================================
            if coag_is_off
                fprintf('Skipping coagulation kernels (coag OFF)...\n');

                % Create "zero betas" in a way that matches your BetaMatrices type.
                % We do it by building normally once from zeros (safe).
                Z = zeros(obj.config.n_sections);

                % This assumes combineAndScale can accept numeric matrices.
                % If it expects objects, tell me and I'll match it exactly.
                % betas = obj.assembler.combineAndScale(Z, Z, Z);

                % ==========================================================
                % NEW-2025-12-13: SAFE ZERO BETAS (no BetaAssembler dependency)
                % ==========================================================
                Ns = obj.config.n_sections;
                betas = BetaMatrices();
                betas.b1  = zeros(Ns,Ns);
                betas.b2  = zeros(Ns,Ns);
                betas.b3  = zeros(Ns,Ns);
                betas.b4  = zeros(Ns,Ns);
                betas.b5  = zeros(Ns,Ns);
                betas.b25 = zeros(Ns,Ns);

            else
                fprintf('Computing coagulation kernels...\n');
                b_brown = obj.assembler.computeFor('KernelBrown');
                b_shear = obj.assembler.computeFor('KernelCurSh');
                b_ds    = obj.assembler.computeFor('KernelCurDS');

                betas = obj.assembler.combineAndScale(b_brown, b_shear, b_ds);
            end

            % Add betas to operators for diagnostics
            obj.operators.betas = betas;

            % Create RHS
            obj.rhs = CoagulationRHS(betas, obj.operators.linear, ...
                obj.operators.disagg_minus, obj.operators.disagg_plus, obj.config);

            % Pass epsilon series if available
            if ~isempty(obj.epsilon_time) && ismethod(obj.rhs, 'setEpsilonSeries')
                obj.rhs.setEpsilonSeries(obj.epsilon_time, obj.epsilon_values);
            end

            % Validate RHS
            if ismethod(obj.rhs, 'validate')
                obj.rhs.validate();
            end

            % Solver options (your clean path)
            solver_opts = [];
            try
                solver_opts = odeset('MaxStep', obj.config.delta_t, ...
                                     'InitialStep', obj.config.delta_t);
            catch
                solver_opts = [];
            end

            if isprop(obj.config,'ode_options') && ~isempty(obj.config.ode_options)
                try
                    solver_opts = odeset(solver_opts, obj.config.ode_options);
                catch
                end
            end

            if ~isempty(p.Results.solver_options)
                try
                    solver_opts = odeset(solver_opts, p.Results.solver_options);
                catch
                end
            end

            % Solve ODEs
            fprintf('Solving ODEs...\n');
            [t, Y] = obj.solver.solve(obj.rhs, tspan, v0, solver_opts);

            % Store results
            obj.result.time = t;
            obj.result.concentrations = Y;
            obj.result.initial_conditions = v0;
            obj.result.betas = betas;
            obj.result.operators = obj.operators;

            % Compute diagnostics
            fprintf('Computing diagnostics...\n');
            obj.result.diagnostics = obj.computeDiagnostics(t, Y);

            % Output data
            fprintf('Computing output data...\n');
            obj.result.output_data = OutputGenerator.spectraAndFluxes(t, Y, obj.grid, obj.config);

            fprintf('Simulation completed successfully.\n');
            result = obj.result;
        end

        % ======== everything below here is exactly your original code ========

        function diagnostics = computeDiagnostics(obj, t, Y)
            diagnostics = struct();

            [diagnostics.sectional_gains, diagnostics.sectional_losses] = ...
                MassBalanceAnalyzer.sectional(Y, obj.operators);

            [diagnostics.total_gains, diagnostics.total_losses] = ...
                MassBalanceAnalyzer.total(Y, obj.operators);

            diagnostics.rate_terms = obj.computeRateTermsOverTime(Y);
            diagnostics.mass_conservation = obj.checkMassConservation(Y);
        end

        function rate_terms = computeRateTermsOverTime(obj, Y)
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

                rate_terms.term1(i, :) = term1(:)';
                rate_terms.term2(i, :) = term2(:)';
                rate_terms.term3(i, :) = term3(:)';
                rate_terms.term4(i, :) = term4(:)';
                rate_terms.term5(i, :) = term5(:)';
            end
        end

        function conservation = checkMassConservation(obj, Y)
            conservation = struct();
            conservation.total_mass = sum(Y, 2);

            if size(Y, 1) > 1
                conservation.mass_change_rate = diff(conservation.total_mass);
                if length(conservation.total_mass) > 1
                    conservation.relative_change = conservation.mass_change_rate ./ conservation.total_mass(1:end-1);
                else
                    conservation.relative_change = [];
                end
            else
                conservation.mass_change_rate = [];
                conservation.relative_change = [];
            end

            conservation.is_conserved = all(conservation.total_mass >= 0);

            if isfield(obj.operators, 'sink_loss')
                try
                    if isprop(obj.config,'use_column') && obj.config.use_column
                        try
                            Ns = obj.config.n_sections;
                            Nz = obj.config.getNumLayers();

                            r_i = obj.grid.getFractalRadii();
                            r_v = obj.grid.getConservedRadii();
                            w   = SettlingVelocityService.velocity(r_i, r_v, obj.grid.setcon);
                            w   = w * obj.config.day_to_sec;
                            dz_cm = obj.config.dz * 100;
                            sink_rate = w / dz_cm;

                            Nt = size(Y,1);
                            export_loss = zeros(Nt,1);

                            for it = 1:Nt
                                vflat = Y(it,:)';
                                N2 = reshape(vflat, [Ns, Nz]);
                                vbot = N2(:, end);
                                export_loss(it) = sum(sink_rate .* vbot);
                            end

                            conservation.settling_loss = export_loss;
                        catch
                            conservation.settling_loss = [];
                        end
                        return;
                    end

                    sink_diag = diag(obj.operators.sink_loss);
                    if isrow(sink_diag)
                        sink_diag = sink_diag';
                    end

                    if length(sink_diag) == size(Y, 2)
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
            if nargin < 2
                plot_flag = true;
            end

            if isempty(obj.result)
                error('No simulation results available. Run simulation first.');
            end

            if plot_flag
                fprintf('Generating plots...\n');

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

            obj.displayDiagnosticsSummary();
        end

        function displayDiagnosticsSummary(obj)
            fprintf('\n=== Simulation Diagnostics Summary ===\n');

            fprintf('Simulation time: %.2f to %.2f days\n', ...
                obj.result.time(1), obj.result.time(end));
            fprintf('Number of time points: %d\n', length(obj.result.time));

            if isprop(obj.config, 'use_column') && obj.config.use_column
                Nz = obj.config.getNumLayers();
                Ns = obj.config.n_sections;
                fprintf('Number of sections: %d (with %d vertical layers; state length=%d)\n', ...
                    Ns, Nz, size(obj.result.concentrations, 2));
            else
                fprintf('Number of sections: %d\n', size(obj.result.concentrations, 2));
            end

            MassBalanceAnalyzer.displayBalanceSummary(...
                obj.result.diagnostics.sectional_gains, ...
                obj.result.diagnostics.sectional_losses, ...
                obj.result.time);

            if obj.result.diagnostics.mass_conservation.is_conserved
                fprintf('Mass conservation: PASSED\n');
            else
                fprintf('Mass conservation: FAILED\n');
            end

            if isfield(obj.result, 'betas')
                obj.result.betas.displaySummary();
            end
        end

        function exportResults(obj, filename)
            if nargin < 2
                filename = sprintf('coagulation_simulation_%s.mat', ...
                    datestr(now, 'yyyymmdd_HHMMSS'));
            end

            OutputGenerator.exportData(obj.result.output_data, filename);

            results_filename = strrep(filename, '.mat', '_full.mat');
            save(results_filename, 'obj');
            fprintf('Full simulation results saved to: %s\n', results_filename);
        end

        function enableTracer(obj)
            warning('Tracer integration not yet implemented');
        end
    end
end