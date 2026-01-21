classdef CoagulationSimulation < handle
    % COAGULATIONSIMULATION Main simulation controller for coagulation system
    %
    % IMPORTANT:
    % - Do NOT build betas or RHS in initializeComponents().
    % - Betas + RHS are built inside run() so cfg switches are respected.
    %
    % This file is intentionally ASCII-only to avoid MATLAB "Invalid text character"
    % parse errors caused by hidden unicode/invisible bytes.
    %
    % NOTE:
    % - This version keeps the same public API you were using:
    %     sim = CoagulationSimulation(cfg); out = sim.run();
    % - It keeps both rateTerms() and decomposeTerms() diagnostics paths.

    properties
        config      % SimulationConfig
        grid        % DerivedGrid
        assembler   % BetaAssembler
        operators   % struct of operators
        rhs         % CoagulationRHS
        solver      % ODESolver
        result      % struct

        % Optional epsilon forcing
        epsilon_time     = [];   % [d]
        epsilon_values   = [];   % scalar eps(t), Nt x 1
        epsilon_values_z = [];   % eps(t,z), Nt x Nz (or Nz x Nt)

        % ==========================================================
        % NEW-2026-01-15: debug helpers (optional)
        % ==========================================================
        debug_epsilon_shapes = false;  % set true if you want prints
    end

    methods
        function obj = CoagulationSimulation(varargin)
            % Constructor: accepts SimulationConfig or parameter-value pairs
            if nargin == 0
                obj.config = SimulationConfig();
            elseif nargin == 1 && isa(varargin{1}, 'SimulationConfig')
                obj.config = varargin{1};
            else
                obj.config = SimulationConfig(varargin{:});
            end

            obj.config.validate();

            try
                fprintf('[INFO] Using CoagulationSimulation from: %s\n', which('CoagulationSimulation'));
            catch
            end

            obj.initializeComponents();
        end

        function initializeComponents(obj)
            % Initialize all simulation components (NO RHS build here)

            obj.grid = obj.config.derive();

            % Grid sanity check
            try
                if isempty(obj.grid)
                    error('Grid is empty after config.derive().');
                end
                if ~(isobject(obj.grid) || isa(obj.grid, 'handle'))
                    error('Grid is not an object (class=%s).', class(obj.grid));
                end
            catch ME
                warning('[WARN] Grid sanity check failed: %s', ME.message);
            end

            obj.assembler = BetaAssembler(obj.config, obj.grid);

            obj.operators = struct();
            obj.operators.config = obj.config;
            obj.operators.grid = obj.grid;

            obj.operators.growth = LinearProcessBuilder.growthMatrix(obj.config, obj.grid);

            [obj.operators.sink_loss, obj.operators.sink_rate, obj.operators.settling_vel_cmday] = ...
                LinearProcessBuilder.sinkingMatrix(obj.config, obj.grid);

            [obj.operators.disagg_minus, obj.operators.disagg_plus] = ...
                LinearProcessBuilder.disaggregationMatrices(obj.config);

            % Prevent double counting if new disagg is applied in RHS
            use_new_disagg = false;
            try
                if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg) && logical(obj.config.enable_disagg)
                    if isprop(obj.config,'disagg_apply_in') && ~isempty(obj.config.disagg_apply_in)
                        use_new_disagg = strcmpi(string(obj.config.disagg_apply_in), "rhs");
                    else
                        use_new_disagg = true;
                    end
                end
            catch
                use_new_disagg = false;
            end

            if use_new_disagg
                try
                    obj.operators.disagg_minus = 0 * obj.operators.disagg_minus;
                    obj.operators.disagg_plus  = 0 * obj.operators.disagg_plus;
                catch
                    Ns = obj.config.n_sections;
                    obj.operators.disagg_minus = zeros(Ns, Ns);
                    obj.operators.disagg_plus  = zeros(Ns, Ns);
                end
            end

            obj.operators.linear = LinearProcessBuilder.linearMatrix(obj.config, obj.grid);

            obj.solver = ODESolver();

            obj.rhs = [];
            obj.result = struct();
        end

        function verifyNoShadowing(obj) %#ok<MANU>
            try
                fprintf('[INFO] which CoagulationSimulation: %s\n', which('CoagulationSimulation'));
                fprintf('[INFO] which SimulationConfig:     %s\n', which('SimulationConfig'));
                fprintf('[INFO] which DerivedGrid:          %s\n', which('DerivedGrid'));
                fprintf('[INFO] which CoagulationRHS:       %s\n', which('CoagulationRHS'));
                fprintf('[INFO] which BetaAssembler:        %s\n', which('BetaAssembler'));
            catch
            end
        end

        function setEpsilonTimeSeries(obj, t, eps_vals)
            obj.epsilon_time = t(:);

            if isvector(eps_vals)
                obj.epsilon_values   = eps_vals(:);
                obj.epsilon_values_z = [];
                if ~isempty(obj.rhs) && ismethod(obj.rhs, 'setEpsilonSeries')
                    obj.rhs.setEpsilonSeries(obj.epsilon_time, obj.epsilon_values);
                end
            else
                obj.epsilon_values   = [];
                obj.epsilon_values_z = eps_vals;
            end
        end

        function setEpsilonTimeSeriesZ(obj, t, eps_mat)
            obj.epsilon_time     = t(:);
            obj.epsilon_values_z = eps_mat;
        end

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
            fprintf('Starting coagulation simulation...\n');

            % Record identity info (keep simple ASCII strings)
            try
                obj.result = struct();
                obj.result.which_CoagulationSimulation = which('CoagulationSimulation');
                obj.result.which_CoagulationRHS        = which('CoagulationRHS');
                obj.result.grid_class = class(obj.grid);
                obj.result.grid_has_av_vol = isprop(obj.grid,'av_vol');
            catch
            end

            % Parse optional args
            p = inputParser;
            addParameter(p, 'tspan', [], @isnumeric);
            addParameter(p, 'v0', [], @isnumeric);
            addParameter(p, 'solver_options', [], @isstruct);
            parse(p, varargin{:});

            % Time span
            if isempty(p.Results.tspan)
                tspan = obj.config.t_init:obj.config.delta_t:obj.config.t_final;
            else
                tspan = p.Results.tspan;
            end

            % Initial condition
            if isempty(p.Results.v0)
                v0 = obj.buildInitialCondition();
            else
                v0 = p.Results.v0;
            end

            % ==========================================================
            % NEW-2026-01-15: normalize cfg epsilon_series shape early
            % - Accept Nz x Nt or Nt x Nz
            % - Convert to Nt x Nz to satisfy SimulationConfig.validate()
            % - Does NOT change physics, only fixes orientation
            % ==========================================================
            try
                if isprop(obj.config,'use_column') && obj.config.use_column
                    if isprop(obj.config,'epsilon_time') && isprop(obj.config,'epsilon_series')
                        if ~isempty(obj.config.epsilon_time) && ~isempty(obj.config.epsilon_series)
                            if ~isvector(obj.config.epsilon_series)
                                Nt = numel(obj.config.epsilon_time(:));
                                Nz = obj.config.getNumLayers();
                                E  = obj.config.epsilon_series;

                                % If E is Nz x Nt, transpose
                                if size(E,1) == Nz && size(E,2) == Nt
                                    if obj.debug_epsilon_shapes
                                        fprintf('[EPS] cfg.epsilon_series is Nz x Nt (%dx%d). Transposing -> Nt x Nz.\n', size(E,1), size(E,2));
                                    end
                                    obj.config.epsilon_series = E.';
                                end

                                % If it is still not Nt x Nz, just print debug
                                if obj.debug_epsilon_shapes
                                    EE = obj.config.epsilon_series;
                                    fprintf('[EPS] cfg.epsilon_series now size=%dx%d (expect Nt=%d, Nz=%d)\n', size(EE,1), size(EE,2), Nt, Nz);
                                end
                            end
                        end
                    end
                end
            catch
            end

            % ==========================================================
            % NEW-2026-01-15: if user uses epsilon_const = NaN, treat as "unset"
            % (some codes interpret NaN as a value; we want it to not override series)
            % ==========================================================
            try
                if isprop(obj.config,'epsilon_const') && ~isempty(obj.config.epsilon_const)
                    if isnumeric(obj.config.epsilon_const) && numel(obj.config.epsilon_const) == 1
                        if isfinite(obj.config.epsilon_const) == 0
                            % leave it as-is, but this is a hint for your RHS:
                            % if RHS checks "isfinite", NaN won't override
                        end
                    end
                end
            catch
            end

            % Mirror cfg epsilon forcing if present and sim has none
            if (isempty(obj.epsilon_time) || (isempty(obj.epsilon_values) && isempty(obj.epsilon_values_z)))
                if isprop(obj.config,'epsilon_time') && isprop(obj.config,'epsilon_series')
                    if ~isempty(obj.config.epsilon_time) && ~isempty(obj.config.epsilon_series)
                        try
                            obj.epsilon_time = obj.config.epsilon_time(:);
                            if ~isvector(obj.config.epsilon_series)
                                obj.epsilon_values_z = obj.config.epsilon_series;
                            else
                                obj.epsilon_values = obj.config.epsilon_series(:);
                            end
                        catch
                        end
                    end
                end
            end

            % If epsilon_values_z exists, push to cfg.epsilon_series (Nt x Nz)
            try
                if ~isempty(obj.epsilon_time) && ~isempty(obj.epsilon_values_z)
                    if isprop(obj.config,'use_column') && obj.config.use_column
                        Nz = obj.config.getNumLayers();
                        E  = obj.epsilon_values_z;

                        if size(E,1) == Nz && size(E,2) == numel(obj.epsilon_time)
                            E = E.'; % Nz x Nt -> Nt x Nz
                        end

                        if size(E,1) == numel(obj.epsilon_time) && size(E,2) == Nz
                            obj.config.epsilon_time   = obj.epsilon_time;
                            obj.config.epsilon_series = E;
                        end
                    end
                end
            catch
            end

            % ==========================================================
            % NEW-2026-01-15: re-run validate after we may have transposed epsilon
            % (keeps your original behavior, but avoids the assert crash)
            % ==========================================================
            try
                obj.config.validate();
            catch ME
                warning('[WARN] config.validate() failed after epsilon normalization: %s', ME.message);
                rethrow(ME);
            end

            % Decide if coagulation is OFF
            coag_is_off = false;
            why_off = {};

            if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag) && ~logical(obj.config.enable_coag)
                coag_is_off = true; why_off{end+1} = 'enable_coag=false';
            end
            if isprop(obj.config,'coag_scale') && ~isempty(obj.config.coag_scale) && obj.config.coag_scale == 0
                coag_is_off = true; why_off{end+1} = 'coag_scale=0';
            end
            if isprop(obj.config,'beta_scale') && ~isempty(obj.config.beta_scale) && obj.config.beta_scale == 0
                coag_is_off = true; why_off{end+1} = 'beta_scale=0';
            end
            if isprop(obj.config,'kernel_scale') && ~isempty(obj.config.kernel_scale) && obj.config.kernel_scale == 0
                coag_is_off = true; why_off{end+1} = 'kernel_scale=0';
            end

            % Disagg status (informational)
            disagg_on = false;
            try
                if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                    disagg_on = logical(obj.config.enable_disagg);
                end
            catch
                disagg_on = false;
            end

            % Column status
            col_on = false;
            Nz_here = 1;
            if isprop(obj.config,'use_column') && obj.config.use_column
                col_on = true;
                Nz_here = obj.config.getNumLayers();
            end

            gm = '';
            try, gm = string(obj.config.growth_mode); catch, gm = ''; end

            if coag_is_off
                why_txt = 'unknown';
                if ~isempty(why_off), why_txt = strjoin(why_off, ', '); end
                fprintf('Physics: coag=OFF (%s) | disagg=%s | growth_mode=%s | column=%s (Nz=%d)\n', ...
                    why_txt, string(disagg_on), gm, string(col_on), Nz_here);
            else
                fprintf('Physics: coag=ON | disagg=%s | growth_mode=%s | column=%s (Nz=%d)\n', ...
                    string(disagg_on), gm, string(col_on), Nz_here);
            end

            % Compute betas
            if coag_is_off
                fprintf('Skipping coagulation kernels (coag OFF)...\n');
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
                betas   = obj.assembler.combineAndScale(b_brown, b_shear, b_ds);
            end

            obj.operators.betas = betas;

            % Create RHS (pass grid if supported)
            try
                obj.rhs = CoagulationRHS(betas, obj.operators.linear, ...
                    obj.operators.disagg_minus, obj.operators.disagg_plus, obj.config, obj.grid);
            catch
                obj.rhs = CoagulationRHS(betas, obj.operators.linear, ...
                    obj.operators.disagg_minus, obj.operators.disagg_plus, obj.config);
                try
                    if ismethod(obj.rhs,'setGrid')
                        obj.rhs.setGrid(obj.grid);
                    end
                catch
                end
            end

            % Pass scalar eps(t) if available
            if ~isempty(obj.epsilon_time) && ~isempty(obj.epsilon_values)
                try
                    if ismethod(obj.rhs,'setEpsilonSeries')
                        obj.rhs.setEpsilonSeries(obj.epsilon_time, obj.epsilon_values);
                    end
                catch
                end
            end

            % Validate RHS
            if ~isempty(obj.rhs) && ismethod(obj.rhs,'validate')
                obj.rhs.validate();
            end

            % Solver options
            solver_opts = [];
            try
                solver_opts = odeset('MaxStep', obj.config.delta_t, 'InitialStep', obj.config.delta_t);
            catch
                solver_opts = [];
            end

            if isprop(obj.config,'ode_options') && ~isempty(obj.config.ode_options)
                try, solver_opts = odeset(solver_opts, obj.config.ode_options); catch, end
            end
            if ~isempty(p.Results.solver_options)
                try, solver_opts = odeset(solver_opts, p.Results.solver_options); catch, end
            end

            try
                if isprop(obj.config,'use_nonnegative') && ~isempty(obj.config.use_nonnegative) && logical(obj.config.use_nonnegative)
                    solver_opts = odeset(solver_opts, 'NonNegative', 1:numel(v0));
                end
            catch
            end

            % Solve
            fprintf('Solving ODEs...\n');
            [t, Y] = obj.solver.solve(obj.rhs, tspan, v0, solver_opts);
            
            if isempty(t) || isempty(Y) || size(Y,1) < 2
                error('ODE solve failed or returned too few time points. Check RHS stability / tolerances.');
            end
            % Store results
            obj.result.time = t;
            obj.result.concentrations = Y;
            obj.result.initial_conditions = v0;
            obj.result.betas = betas;
            obj.result.operators = obj.operators;

            obj.result.flags = struct();
            obj.result.flags.coag_is_off = coag_is_off;
            obj.result.flags.disagg_on   = disagg_on;
            obj.result.flags.column_on   = col_on;

            % Diagnostics
            fprintf('Computing diagnostics...\n');
            obj.result.diagnostics = obj.computeDiagnostics(t, Y);

            % Output data
            fprintf('Computing output data...\n');
            obj.result.output_data = OutputGenerator.spectraAndFluxes(t, Y, obj.grid, obj.config);

            fprintf('Simulation completed successfully.\n');
            result = obj.result;
        end

        % ===== diagnostics =====

        function diagnostics = computeDiagnostics(obj, t, Y) %#ok<INUSD>
    diagnostics = struct();

    % ------------------------------------------------------------
    % NEW: robust MassBalanceAnalyzer calls (do NOT crash run)
    % ------------------------------------------------------------
    try
        % Path sanity (helps catch shadowing / wrong file)
        try
            diagnostics.which_MassBalanceAnalyzer = which('MassBalanceAnalyzer');
        catch
            diagnostics.which_MassBalanceAnalyzer = '';
        end

        % Check class + method existence
        hasClass = (exist('MassBalanceAnalyzer','class') == 8);
        hasSect  = false;
        hasTotal = false;

        if hasClass
            try
                m = methods('MassBalanceAnalyzer');
                hasSect  = any(strcmp(m,'sectional'));
                hasTotal = any(strcmp(m,'total'));
            catch
                hasSect  = false;
                hasTotal = false;
            end
        end

        % ---- sectional ----
        if hasClass && hasSect
            [diagnostics.sectional_gains, diagnostics.sectional_losses] = ...
                MassBalanceAnalyzer.sectional(Y, obj.operators);
        else
            % OLD line (kept): MassBalanceAnalyzer.sectional(Y, obj.operators);
            warning('computeDiagnostics: MassBalanceAnalyzer.sectional not found. Skipping sectional diagnostics.');
            diagnostics.sectional_gains  = struct();
            diagnostics.sectional_losses = struct();
        end

        % ---- total ----
        if hasClass && hasTotal
            [diagnostics.total_gains, diagnostics.total_losses] = ...
                MassBalanceAnalyzer.total(Y, obj.operators);
        else
            % OLD line (kept): MassBalanceAnalyzer.total(Y, obj.operators);
            warning('computeDiagnostics: MassBalanceAnalyzer.total not found. Skipping total diagnostics.');
            diagnostics.total_gains  = struct();
            diagnostics.total_losses = struct();
        end

    catch ME
        warning('computeDiagnostics: MassBalanceAnalyzer failed (%s). Continuing without these diagnostics.', ME.message);
        diagnostics.sectional_gains  = struct();
        diagnostics.sectional_losses = struct();
        diagnostics.total_gains      = struct();
        diagnostics.total_losses     = struct();
        diagnostics.massbalance_error = ME.message;
    end

    % Keep the rest exactly as before
    diagnostics.rate_terms = obj.computeRateTermsOverTime(Y);
    diagnostics.mass_conservation = obj.checkMassConservation(Y);
end

        function rate_terms = computeRateTermsOverTime(obj, Y)
            n_times = size(Y, 1);
            n_state = size(Y, 2);

            rate_terms = struct();
            rate_terms.term1 = zeros(n_times, n_state);
            rate_terms.term2 = zeros(n_times, n_state);
            rate_terms.term3 = zeros(n_times, n_state);
            rate_terms.term4 = zeros(n_times, n_state);
            rate_terms.term5 = zeros(n_times, n_state);

            rate_terms.dv_lin    = NaN(n_times, n_state);
            rate_terms.dv_pp     = NaN(n_times, n_state);
            rate_terms.dv_disagg = NaN(n_times, n_state);
            rate_terms.dv_coag   = NaN(n_times, n_state);
            rate_terms.dv_tot    = NaN(n_times, n_state);

            for i = 1:n_times
                v = Y(i, :)';

                % legacy rateTerms if present
                try
                    if ~isempty(obj.rhs) && ismethod(obj.rhs,'rateTerms')
                        [term1, term2, term3, term4, term5] = obj.rhs.rateTerms(v);
                        rate_terms.term1(i, :) = term1(:)';
                        rate_terms.term2(i, :) = term2(:)';
                        rate_terms.term3(i, :) = term3(:)';
                        rate_terms.term4(i, :) = term4(:)';
                        rate_terms.term5(i, :) = term5(:)';
                    end
                catch
                end

                % decomposeTerms if present
                try
                    if ~isempty(obj.rhs) && ismethod(obj.rhs,'decomposeTerms')
                        tt = 0;
                        if isfield(obj.result,'time') && numel(obj.result.time) >= i
                            tt = obj.result.time(i);
                        end
                        T = obj.rhs.decomposeTerms(tt, v);

                        rate_terms.dv_tot(i,:)    = T.dv_tot(:)';
                        rate_terms.dv_lin(i,:)    = T.dv_lin(:)';
                        rate_terms.dv_pp(i,:)     = T.dv_pp(:)';
                        rate_terms.dv_disagg(i,:) = T.dv_disagg(:)';
                        rate_terms.dv_coag(i,:)   = T.dv_coag(:)';
                    end
                catch
                end
            end
        end

        function conservation = checkMassConservation(obj, Y)
            conservation = struct();
            conservation.total_mass = sum(Y, 2);

            if size(Y, 1) > 1
                conservation.mass_change_rate = diff(conservation.total_mass);
                if numel(conservation.total_mass) > 1
                    conservation.relative_change = conservation.mass_change_rate ./ conservation.total_mass(1:end-1);
                else
                    conservation.relative_change = [];
                end
            else
                conservation.mass_change_rate = [];
                conservation.relative_change = [];
            end

            conservation.is_conserved = all(conservation.total_mass >= 0);

            conservation.settling_loss = [];
            if isfield(obj.operators, 'sink_loss')
                try
                    sink_diag = diag(obj.operators.sink_loss);
                    if isrow(sink_diag), sink_diag = sink_diag'; end
                    if numel(sink_diag) == size(Y, 2)
                        conservation.settling_loss = Y * sink_diag;
                    end
                catch
                end
            end
        end

        function exportResults(obj, filename)
            if nargin < 2
                filename = sprintf('coagulation_simulation_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
            end

            OutputGenerator.exportData(obj.result.output_data, filename);

            results_filename = strrep(filename, '.mat', '_full.mat');
            save(results_filename, 'obj');
            fprintf('Full simulation results saved to: %s\n', results_filename);
        end
    end
end
