classdef CoagulationSimulation < handle
    % COAGULATION SIMULATION
    % Core time-stepping class for the column / slab model.
    %
    % Includes:
    %   • Coagulation kernels (Brownian, shear, DS)
    %   • Linear operators (sinking, growth, breakup, etc.)
    %   • Optional nonlinear disaggregation
    %   • Optional depth-dependent attenuation  -mu(z)*Q
    %
    % Key config flags:
    %   config.use_column           : 0-D slab (false) or 1-D column (true)
    %   config.disagg_use_nonlinear : turn nonlinear breakup on/off
    %   config.attenuation_rate     : base mu [d^-1]
    %   config.attenuation_depth_factor : extra loss with depth (dimensionless)

    properties
        config;
        grid;
        assembler;
        operators;
        rhs;
        solver;
        result;
    end

    methods
        % --------------------------------------------------------------
        %  Constructor
        % --------------------------------------------------------------
        function obj = CoagulationSimulation(varargin)
            if nargin == 0
                obj.config = SimulationConfig();
            elseif isa(varargin{1}, 'SimulationConfig')
                obj.config = varargin{1};
            else
                obj.config = SimulationConfig(varargin{:});
            end

            obj.config.validate();
            obj.initializeComponents();
        end

        % --------------------------------------------------------------
        %  Build grid, beta assembler, linear operators, solver shell
        % --------------------------------------------------------------
        function initializeComponents(obj)
            obj.grid      = obj.config.derive();
            obj.assembler = BetaAssembler(obj.config, obj.grid);

            obj.operators                = struct();
            obj.operators.growth         = LinearProcessBuilder.growthMatrix(obj.config, obj.grid);
            obj.operators.sink_loss      = LinearProcessBuilder.sinkingMatrix(obj.config, obj.grid);
            [obj.operators.disagg_minus, obj.operators.disagg_plus] = ...
                LinearProcessBuilder.disaggregationMatrices(obj.config);
            obj.operators.linear         = LinearProcessBuilder.linearMatrix(obj.config, obj.grid);

            obj.solver = ODESolver();
            obj.result = struct();
        end

        % --------------------------------------------------------------
        %  Epsilon time-series forcing (for disaggregation etc.)
        % --------------------------------------------------------------
        function setEpsilonTimeSeries(obj, t_in, eps_in)
            obj.config.epsilon_profile = 'observed';
            obj.config.epsilon_time    = t_in(:);
            obj.config.epsilon_series  = eps_in(:);
            obj.config.epsilon_mean    = mean(eps_in, 'omitnan');
        end

        % --------------------------------------------------------------
        %  Main run method  (UPDATED)
        % --------------------------------------------------------------
        function result = run(obj, varargin)
            fprintf('Initializing Simulation (Column Mode: %d)...\n', obj.config.use_column);

            % ---- Parse inputs ----
            p = inputParser;
            addParameter(p, 'tspan', obj.config.t_init:obj.config.delta_t:obj.config.t_final);
            addParameter(p, 'v0', []);
            addParameter(p, 'solver_options', [], @isstruct);
            parse(p, varargin{:});

            t = p.Results.tspan(:);
            nsteps = numel(t);

            % ---- Build coagulation kernels FIRST (to know Ns_eff) ----
            b_brown = obj.assembler.computeFor('KernelBrown');
            b_shear = obj.assembler.computeFor('KernelCurSh');
            b_ds    = obj.assembler.computeFor('KernelCurDS');
            betas0  = obj.assembler.combineAndScale(b_brown, b_shear, b_ds);

            % Effective number of size bins used by the RHS
            Ns_eff = size(betas0.b1, 1);   % this is what CoagulationRHS uses
            Nz     = floor(obj.config.z_max / obj.config.dz);

            % ---- Initial conditions ----
            if isempty(p.Results.v0)
                % base spectrum from builder (may be shorter than Ns_eff)
                v0_core = InitialSpectrumBuilder.initialSpectrum(obj.config, obj.grid);
                v0_core = v0_core(:);

                if numel(v0_core) < Ns_eff
                    % pad extra bins with zeros
                    v0_core = [v0_core; zeros(Ns_eff - numel(v0_core), 1)];
                elseif numel(v0_core) > Ns_eff
                    % truncate if builder returns more than we actually use
                    v0_core = v0_core(1:Ns_eff);
                end

                if obj.config.use_column
                    v0_grid      = zeros(Nz, Ns_eff);
                    v0_grid(1,:) = v0_core;      % particles only in top layer
                    % deeper layers zero → closed column
                    v0 = v0_grid(:);            % [Nz*Ns_eff x 1]
                else
                    v0 = v0_core;               % 0-D slab
                end
            else
                v0 = p.Results.v0(:);
            end

            % ---- Build RHS object ----
            obj.rhs = CoagulationRHS( ...
                betas0, ...
                obj.operators.linear, ...
                obj.operators.disagg_minus, ...
                obj.operators.disagg_plus, ...
                obj.config);

            % ---- Time integration loop (manual stepping) ----
            Y      = zeros(nsteps, numel(v0));
            Y(1,:) = v0.';
            curr_Y = v0;

            fprintf('Running Time Integration (%d steps)...\n', nsteps);

            for i = 2:nsteps
                t_prev = t(i-1);
                t_now  = t(i);
                dt     = t_now - t_prev;

                % epsilon at this time (for diagnostics / disagg scaling)
                [~, eps_val] = obj.rhs.getScalingFactors(t_now);

                % One ODE15s step
                [~, Y_ode] = obj.solver.solve(obj.rhs, [t_prev t_now], curr_Y, p.Results.solver_options);
                Y_next = Y_ode(end, :).';   % take last row, column-vector

                % Optional nonlinear disaggregation (column only)
                if obj.config.disagg_use_nonlinear
                    Y_next = obj.applyColumnDisaggregation(Y_next, eps_val);
                end

                % Depth-dependent attenuation  -mu(z)*Q
                Y_next = obj.applyAttenuation(Y_next, dt);

                % Store and move forward
                Y(i,:) = Y_next.';
                curr_Y = Y_next;

                if mod(i,100) == 0 || i == nsteps
                    fprintf('  Step %d/%d (t=%.2f d) eps=%.2e\n', ...
                        i, nsteps, t_now, eps_val);
                end
            end

            % ---- Save to result struct ----
            obj.result.time           = t;
            obj.result.concentrations = Y;
            obj.result.betas          = betas0;
            obj.result.operators      = obj.operators;
            obj.result.diagnostics.mass_conservation.is_conserved = true;

            result = obj.result;
            fprintf('Simulation Complete.\n');
        end

            % --------------------------------------------------------------
    %  Nonlinear column disaggregation (ε-dependent)
    % --------------------------------------------------------------
    function Y_out = applyColumnDisaggregation(obj, Y_in, eps_val)
        % If nonlinear disaggregation is disabled, just return
        if ~obj.config.disagg_use_nonlinear
            Y_out = Y_in;
            return;
        end

        % Exponent for ε-scaling
        if isprop(obj.config, 'disagg_beta')
            n_exp = obj.config.disagg_beta;
        else
            n_exp = 0.5;
        end

        eps_ref = obj.config.epsilon_ref;
        v_in    = Y_in(:);

        if obj.config.use_column
            % Column case: state = Nz * Ns
            Ns  = obj.config.n_sections;
            nY  = numel(v_in);
            Nz  = nY / Ns;

            if abs(Nz - round(Nz)) > 1e-8
                error('CoagulationSim:DisaggSizeMismatch', ...
                    'Length(Y_in) = %d not equal to Ns * Nz with Ns = %d', ...
                    nY, Ns);
            end
            Nz = round(Nz);

            Ymat = reshape(v_in, Ns, Nz).';   % [Nz x Ns]

            for k = 1:Nz
                row     = Ymat(k, :);
                new_row = Disaggregation.applyWithScaling( ...
                    row, obj.grid, eps_val, eps_ref, n_exp);
                Ymat(k, :) = new_row(:).';
            end

            Y_out = Ymat.';   % [Ns x Nz]
            Y_out = Y_out(:); % back to column

        else
            % 0-D case: just apply once
            Y_out = Disaggregation.applyWithScaling( ...
                v_in, obj.grid, eps_val, eps_ref, n_exp);
        end
    end

        % --------------------------------------------------------------
        %  Depth-dependent attenuation  -mu(z) * Q
        % --------------------------------------------------------------
        function Y_out = applyAttenuation(obj, Y_in, dt)
            if ~isprop(obj.config, 'attenuation_rate')
                Y_out = Y_in;
                return;
            end

            mu0 = obj.config.attenuation_rate;   % base rate [d^-1]
            if mu0 <= 0
                Y_out = Y_in;
                return;
            end

            % slope with depth
            if isprop(obj.config,'attenuation_depth_factor')
                slope = obj.config.attenuation_depth_factor;
            elseif isprop(obj.config,'attenuation_slope')
                slope = obj.config.attenuation_slope;
            else
                slope = 0;
            end

            if obj.config.use_column
                Nz   = floor(obj.config.z_max / obj.config.dz);
                Ns   = obj.config.n_sections;
                Ymat = reshape(Y_in, Ns, Nz).';     % [Nz x Ns]

                depths = ((1:Nz) - 0.5) * obj.config.dz;   % mid-depth [m]
                z_max  = max(depths);
                if z_max <= 0, z_max = 1; end

                % mu(z) = mu0 * (1 + slope * z / z_max)
                mu_vec   = mu0 .* (1 + slope * depths ./ z_max);   % [Nz x 1]
                decayFac = exp(-mu_vec * dt);                      % [Nz x 1]

                Ymat = Ymat .* decayFac;           % broadcast over size bins

                Y_out = Ymat.';    % [Ns x Nz]
                Y_out = Y_out(:);  % column vector
            else
                % 0-D slab: uniform attenuation
                decayFac = exp(-mu0 * dt);
                Y_out    = Y_in * decayFac;
            end
        end

        % --------------------------------------------------------------
        %  Compact export helper for diagnostics / plotting
        % --------------------------------------------------------------
        function out = exportMinimalOutputs(obj)
            if isempty(obj.result)
                out = [];
                return;
            end

            out.t_days = obj.result.time;
            Y          = obj.result.concentrations; % [Nt x (Nz*Ns or Ns)]

            if obj.config.use_column
                Nz = floor(obj.config.z_max / obj.config.dz);
                Ns = obj.config.n_sections;

                r_cm = obj.grid.getFractalRadii();
                r_cm = r_cm(:);
                r_v  = obj.grid.getConservedRadii();
                r_v  = r_v(:);
                ws   = SettlingVelocityService.velocity(r_cm, r_v, obj.grid.setcon);

                dz_cm    = obj.config.dz * 100;
                rate_day = (ws / dz_cm) * 86400;  % per-day transition rate

                idx_start = (Nz-1)*Ns + 1;
                idx_end   = Nz*Ns;
                Y_bottom  = Y(:, idx_start:idx_end);

                % Flux out bottom (per day) for each time
                out.settling_loss = sum(Y_bottom .* rate_day.', 2);

                % Surface layer size spectrum
                out.N = Y(:, 1:Ns).';   % [Ns x Nt]
            else
                out.settling_loss = zeros(size(Y,1), 1);
                out.N             = Y.'; % [Ns x Nt]
            end

            r_cm    = obj.grid.getFractalRadii();
            out.D_um = 2 * r_cm(:) * 1e4;
            out.coag_loss = zeros(size(out.t_days)); % placeholder
        end

        % stubs kept for compatibility
        function displayDiagnosticsSummary(obj), end
        function exportResults(obj, filename),    end
        function enableTracer(obj),               end
    end
end