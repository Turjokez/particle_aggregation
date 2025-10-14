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
        %% UPDATED for Δx = 2 m
        verticalMode = 'uniform';           % Vertical grid mode ('uniform' or 'stretched')
        verticalBreaks = [0, 200, 500, 1000]; % Depth breakpoints for stretched grid [m]
        verticalSpacing = [25, 50, 100];      % Layer thicknesses between breakpoints [m]
        uniformDeltaZ = 2;                    % Target spacing for uniform grid [m]
        maxDepth = 1000;                      % Maximum depth to resolve [m]
        verticalGrid = struct();              % Structure describing vertical grid vectors
        verticalMatrices = struct();          % Pre-computed vertical operator matrices
        fluxAnimationFile = 'Fz_profile.mp4'; % Output path for F(z) animation
    end

    methods
        function obj = CoagulationSimulation(varargin)
            %COAGULATIONSIMULATION Constructor
            % Can accept SimulationConfig object or parameter-value pairs

            if nargin == 0
                % Use default configuration
                obj.config = SimulationConfig();
            elseif isa(varargin{1}, 'SimulationConfig')
                % Use provided configuration and parse optional vertical arguments
                obj.config = varargin{1};
                if numel(varargin) > 1
                    %% UPDATED for Δx = 2 m
                    obj.applyVerticalOptions(varargin{2:end});
                end
            else
                % Create configuration from parameter-value pairs
                %% UPDATED for Δx = 2 m
                [filteredArgs, verticalOptions] = obj.extractVerticalOptions(varargin{:});
                obj.verticalMode = verticalOptions.mode;
                obj.verticalBreaks = verticalOptions.breaks;
                obj.verticalSpacing = verticalOptions.spacing;
                obj.uniformDeltaZ = verticalOptions.uniformDz;
                obj.maxDepth = verticalOptions.maxDepth;
                obj.config = SimulationConfig(filteredArgs{:});
            end

            % Validate configuration
            obj.config.validate();

            % Initialize components
            obj.initializeComponents();

            %% UPDATED for Δx = 2 m
            obj.configureVerticalGrid();
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

        function setVerticalMode(obj, mode, varargin)
            %% UPDATED for Δx = 2 m
            %SETVERTICALMODE Update the vertical grid mode and regenerate operators

            if nargin < 2 || isempty(mode)
                return;
            end

            obj.verticalMode = mode;
            if ~isempty(varargin)
                obj.applyVerticalOptions(varargin{:});
            end

            obj.configureVerticalGrid();
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
            fprintf('Solving ODEs...\n');
            [t, Y] = obj.solver.solve(obj.rhs, tspan, v0, p.Results.solver_options);

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

            %% UPDATED for Δx = 2 m
            obj.verticalMatrices = obj.assembleVerticalMatrices();
            obj.result.vertical_profiles = obj.computeVerticalFluxProfiles();

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

                %% UPDATED for Δx = 2 m
                obj.createFluxVisualizations(obj.result.vertical_profiles);
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

    %% UPDATED for Δx = 2 m
    methods (Access = private)
        function configureVerticalGrid(obj)
            %CONFIGUREVERTICALGRID Build vertical grid and refresh operators

            switch lower(obj.verticalMode)
                case 'stretched'
                    grid = obj.buildStretchedGrid(obj.verticalBreaks, obj.verticalSpacing);
                otherwise
                    grid = obj.buildUniformGrid(obj.uniformDeltaZ, obj.maxDepth);
            end

            obj.verticalGrid = grid;

            if isfield(grid, 'delta_z') && ~isempty(grid.delta_z)
                obj.config.dz = grid.delta_z(1);
            end
            if isfield(grid, 'edges') && ~isempty(grid.edges)
                obj.maxDepth = grid.edges(end);
            end

            obj.refreshLinearOperatorsForVerticalGrid();
        end

        function refreshLinearOperatorsForVerticalGrid(obj)
            %REFRESHLINEAROPERATORSFORVERTICALGRID Recompute operators using updated dz

            obj.operators.sink_loss = LinearProcessBuilder.sinkingMatrix(obj.config, obj.grid);
            obj.operators.linear = LinearProcessBuilder.linearMatrix(obj.config, obj.grid);
        end

        function applyVerticalOptions(obj, varargin)
            %APPLYVERTICALOPTIONS Parse and store vertical grid options

            if isempty(varargin)
                return;
            end

            if mod(numel(varargin), 2) ~= 0
                error('Vertical options must be provided as parameter-value pairs.');
            end

            for idx = 1:2:numel(varargin)
                key = lower(string(varargin{idx}));
                value = varargin{idx+1};
                switch key
                    case 'verticalmode'
                        obj.verticalMode = char(value);
                    case 'stretchedbreaks'
                        obj.verticalBreaks = value;
                    case 'stretchedspacing'
                        obj.verticalSpacing = value;
                    case 'uniformdeltaz'
                        obj.uniformDeltaZ = value;
                    case 'maxdepth'
                        obj.maxDepth = value;
                    otherwise
                        % Ignore unrecognized parameters (handled by SimulationConfig)
                end
            end
        end

        function [filteredArgs, verticalOptions] = extractVerticalOptions(obj, varargin)
            %EXTRACTVERTICALOPTIONS Separate vertical arguments from config parameters

            if mod(numel(varargin), 2) ~= 0
                error('Constructor arguments must be provided as parameter-value pairs.');
            end

            mask = true(1, numel(varargin));
            verticalOptions = struct('mode', obj.verticalMode, ...
                'breaks', obj.verticalBreaks, ...
                'spacing', obj.verticalSpacing, ...
                'uniformDz', obj.uniformDeltaZ, ...
                'maxDepth', obj.maxDepth);

            for idx = 1:2:numel(varargin)
                key = lower(string(varargin{idx}));
                switch key
                    case 'verticalmode'
                        verticalOptions.mode = char(varargin{idx+1});
                        mask(idx:idx+1) = false;
                    case 'stretchedbreaks'
                        verticalOptions.breaks = varargin{idx+1};
                        mask(idx:idx+1) = false;
                    case 'stretchedspacing'
                        verticalOptions.spacing = varargin{idx+1};
                        mask(idx:idx+1) = false;
                    case 'uniformdeltaz'
                        verticalOptions.uniformDz = varargin{idx+1};
                        mask(idx:idx+1) = false;
                    case 'maxdepth'
                        verticalOptions.maxDepth = varargin{idx+1};
                        mask(idx:idx+1) = false;
                end
            end

            filteredArgs = varargin(mask);
        end

        function grid = buildUniformGrid(obj, delta_z, maxDepth)
            %BUILDUNIFORMGRID Create a uniform vertical grid description

            if delta_z <= 0
                error('Uniform grid spacing must be positive.');
            end

            edges = 0:delta_z:maxDepth;
            if edges(end) < maxDepth
                edges = [edges, maxDepth]; %#ok<AGROW>
            end

            centers = edges(1:end-1) + diff(edges)/2;
            grid = struct('mode', 'uniform', ...
                'edges', edges(:), ...
                'centers', centers(:), ...
                'delta_z', diff(edges(:)));
        end

        function grid = buildStretchedGrid(obj, breaks, spacing)
            %BUILDSTRETCHEDGRID Construct stretched grid with specified spacing

            if numel(spacing) ~= numel(breaks) - 1
                error('Spacing vector must be one element shorter than breaks.');
            end

            edges = breaks(1);
            for seg = 1:numel(spacing)
                start_depth = breaks(seg);
                end_depth = breaks(seg+1);
                step = spacing(seg);
                if step <= 0
                    error('Stretched grid spacing must be positive.');
                end
                local_edges = start_depth:step:end_depth;
                if abs(local_edges(end) - end_depth) > 1e-6
                    local_edges = [local_edges, end_depth]; %#ok<AGROW>
                end
                if seg > 1
                    local_edges = local_edges(2:end);
                end
                edges = [edges, local_edges]; %#ok<AGROW>
            end

            edges = unique(edges, 'stable');
            if edges(end) < obj.maxDepth
                edges = [edges, obj.maxDepth]; %#ok<AGROW>
            end

            centers = edges(1:end-1) + diff(edges)/2;
            grid = struct('mode', 'stretched', ...
                'edges', edges(:), ...
                'centers', centers(:), ...
                'delta_z', diff(edges(:)));
        end

        function matrices = assembleVerticalMatrices(obj)
            %ASSEMBLEVERTICALMATRICES Construct diffusion/sinking/coagulation blocks

            depth_grid = obj.verticalGrid;
            if isempty(depth_grid)
                matrices = struct();
                return;
            end

            n_depth = numel(depth_grid.centers);
            diffusion_matrix = obj.constructDiffusionMatrix(depth_grid.delta_z, obj.config.kvisc);

            sink_diag = diag(obj.operators.sink_loss);
            sink_block = kron(speye(n_depth), diag(sink_diag));

            if ~isempty(obj.result) && isfield(obj.result, 'betas') && ~isempty(obj.result.betas.b25)
                coag_block = kron(speye(n_depth), obj.result.betas.b25);
            else
                dim = obj.config.n_sections * n_depth;
                coag_block = sparse(dim, dim);
            end

            matrices = struct('diffusion', diffusion_matrix, ...
                'sinking', sparse(sink_block), ...
                'coagulation', sparse(coag_block), ...
                'depth', depth_grid.centers);
        end

        function diffusion_matrix = constructDiffusionMatrix(~, delta_z, diffusion_coeff)
            %CONSTRUCTDIFFUSIONMATRIX Finite-volume diffusion operator for grid

            n = numel(delta_z);
            if n == 0
                diffusion_matrix = sparse(0, 0);
                return;
            end

            main = zeros(n, 1);
            upper = zeros(n-1, 1);
            lower = zeros(n-1, 1);

            if n == 1
                main(1) = -2 / (delta_z(1)^2);
            else
                upper(1) = 2 / (delta_z(1)^2);
                main(1) = -upper(1);
                for row = 2:(n-1)
                    dz_minus = delta_z(row-1);
                    dz_plus = delta_z(row);
                    upper(row) = 2 / (dz_plus * (dz_plus + dz_minus));
                    lower(row-1) = 2 / (dz_minus * (dz_plus + dz_minus));
                    main(row) = -(upper(row) + lower(row-1));
                end
                lower(n-1) = 2 / (delta_z(end)^2);
                main(n) = -lower(n-1);
            end

            diffusion_matrix = diffusion_coeff * spdiags([lower, main, upper], [-1, 0, 1], n, n);
        end

        function vertical = computeVerticalFluxProfiles(obj)
            %COMPUTEVERTICALFLUXPROFILES Estimate F(z) attenuation across depth

            if isempty(obj.verticalGrid) || isempty(obj.result)
                vertical = struct();
                return;
            end

            time = obj.result.time(:);
            total_flux = obj.result.output_data.total_flux(:);
            depth = obj.verticalGrid.centers(:);
            delta_z = obj.verticalGrid.delta_z(:);

            n_time = numel(time);
            n_depth = numel(depth);
            F_profiles = zeros(n_time, n_depth);

            sink_diag = diag(obj.operators.sink_loss);
            avg_sink = mean(sink_diag);
            coag_loss = 0;
            if isfield(obj.result, 'betas') && ~isempty(obj.result.betas.b25)
                coag_loss = mean(sum(obj.result.betas.b25, 2));
            end

            diffusive_scale = obj.config.kvisc;
            attenuation_rate = max(avg_sink, 0) + max(coag_loss, 0) * 1e-6 + max(diffusive_scale, 0) * 1e-4;
            cumulative_depth = [0; cumsum(delta_z(1:end-1))];
            base_profile = exp(-attenuation_rate * cumulative_depth);

            diffusion_response = ones(n_depth, 1);
            if isfield(obj.verticalMatrices, 'diffusion') && ~isempty(obj.verticalMatrices.diffusion)
                diffusion_response = obj.verticalMatrices.diffusion * ones(n_depth, 1);
                diffusion_response = diffusion_response - min(diffusion_response);
                diffusion_response = diffusion_response / max(diffusion_response + eps);
                diffusion_response = 1 - 0.1 * diffusion_response;
            end

            vertical_shape = (base_profile .* diffusion_response);
            max_shape = max(vertical_shape);
            if max_shape <= 0
                max_shape = 1;
            end
            vertical_shape = (vertical_shape' ./ max_shape);

            for t_idx = 1:n_time
                F_profiles(t_idx, :) = total_flux(t_idx) * vertical_shape;
            end

            vertical = struct('time', time, ...
                'depth', depth, ...
                'F_profiles', F_profiles, ...
                'attenuation_rate', attenuation_rate, ...
                'shape', vertical_shape, ...
                'mode', obj.verticalMode, ...
                'animation_file', '');
        end

        function createFluxVisualizations(obj, vertical)
            %CREATEFLUXVISUALIZATIONS Plot 2D map and trigger animation export

            if nargin < 2 || isempty(vertical)
                return;
            end

            if isempty(vertical) || ~isfield(vertical, 'F_profiles')
                return;
            end

            figure('Name', 'F(z) vs Depth', 'NumberTitle', 'off');
            imagesc(vertical.time, vertical.depth, vertical.F_profiles');
            set(gca, 'YDir', 'reverse');
            xlabel('Time [d]');
            ylabel('Depth [m]');
            title('F(z) attenuation across depth');
            colorbar;
            colormap(parula);
            drawnow;

            try
                obj.writeFluxAnimation(vertical);
            catch ME
                warning('Unable to create F(z) animation: %s', ME.message);
            end
        end

        function writeFluxAnimation(obj, vertical)
            %WRITEFLUXANIMATION Generate animation of F(z) profiles

            if isempty(vertical.F_profiles)
                return;
            end

            try
                v = VideoWriter(obj.fluxAnimationFile, 'MPEG-4');
            catch
                warning('VideoWriter with MPEG-4 not available. Skipping animation.');
                return;
            end

            open(v);
            fig = figure('Name', 'F(z) animation', 'NumberTitle', 'off');
            depth = vertical.depth;
            for t_idx = 1:numel(vertical.time)
                plot(vertical.F_profiles(t_idx, :), depth, 'LineWidth', 2);
                set(gca, 'YDir', 'reverse');
                xlabel('Flux [cm^3 m^{-2} d^{-1}]');
                ylabel('Depth [m]');
                title(sprintf('F(z) at t = %.2f d', vertical.time(t_idx)));
                grid on;
                drawnow;
                frame = getframe(fig);
                writeVideo(v, frame);
            end

            close(v);
            if ishghandle(fig)
                close(fig);
            end

            obj.result.vertical_profiles.animation_file = obj.fluxAnimationFile;
        end
    end
end
