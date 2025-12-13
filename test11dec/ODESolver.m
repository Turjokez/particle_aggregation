classdef ODESolver < handle
    % ODESOLVER Adapter for MATLAB ODE solvers (robust + column-safe)

    properties
        solver_type = 'ode15s';  % Default solver
        options;                 % Base options (AbsTol built inside solve)
    end

    methods
        function obj = ODESolver(varargin)
            if nargin > 0
                obj.solver_type = varargin{1};
            end
            obj.setDefaultOptions();
        end

        function [t, Y] = solve(obj, rhs, tspan, v0, custom_options)
            %SOLVE Solve ODE system
            % rhs must implement evaluate(t,v). jacobian(t,v) optional.

            % -----------------------------
            % Basic interface checks
            % -----------------------------
            if isempty(rhs)
                error('ODESolver: RHS is empty.');
            end
            if ~ismethod(rhs, 'evaluate')
                error('ODESolver: RHS must implement evaluate(t, v).');
            end
            if ismethod(rhs, 'validate')
                rhs.validate();
            end

            % Force v0 to column
            v0 = v0(:);
            n_state = numel(v0);

            % -----------------------------
            % Build AbsTol vector to match state length
            % -----------------------------
            abs_tol = 1.0e-18;
            calcomp = (1:n_state);
            at = (abs_tol * 1.5 .^ (-(calcomp-1)));
            at = at(:);

            % -----------------------------
            % Base solver options
            % -----------------------------
            options = obj.options;

            % Use requested output times if tspan is a vector grid
            use_output_times = numel(tspan) > 2;

            % Reasonable strict tolerance
            options = odeset(options, ...
                'RelTol', 1.0e-12, ...
                'AbsTol', at, ...
                'Refine', 0);

            if use_output_times
                dt_out = min(diff(tspan));
                if isfinite(dt_out) && dt_out > 0
                    options = odeset(options, 'MaxStep', dt_out);
                end
            end

            % -----------------------------
            % Merge any custom options
            % -----------------------------
            if nargin >= 5 && ~isempty(custom_options)
                options = odeset(options, custom_options);
            end

            % ----------------------------------------------------------
            % OLD (kept): always forcing NonNegative can break budgets
            % ----------------------------------------------------------
            % options = odeset(options, 'NonNegative', 1:n_state);

            % ----------------------------------------------------------
            % NEW: respect cfg.use_nonnegative if available
            % default: true (for normal runs), but tests can turn it off
            % ----------------------------------------------------------
            use_nn = true;
            try
                if isprop(rhs,'config') && isprop(rhs.config,'use_nonnegative')
                    use_nn = logical(rhs.config.use_nonnegative);
                end
            catch
                use_nn = true;
            end

            if use_nn
                options = odeset(options, 'NonNegative', 1:n_state);
            else
                % OLD (commented): options = odeset(options, 'NonNegative', []);
                % NEW: do nothing (leave NonNegative unset)
            end

            % -----------------------------
            % Jacobian hook (optional)
            % -----------------------------
            if ismethod(rhs, 'jacobian')
                options = odeset(options, 'Jacobian', @(t,v) rhs.jacobian(t,v));
            end

            % RHS wrapper
            f_rhs = @(t,v) rhs.evaluate(t,v);

            % -----------------------------
            % Solve
            % -----------------------------
            switch obj.solver_type
                case 'ode15s'
                    if use_output_times
                        [t, Y] = ode15s(f_rhs, tspan(:), v0, options);
                    else
                        [t, Y] = ode15s(f_rhs, tspan, v0, options);
                    end
                case 'ode23s'
                    if use_output_times
                        [t, Y] = ode23s(f_rhs, tspan(:), v0, options);
                    else
                        [t, Y] = ode23s(f_rhs, tspan, v0, options);
                    end
                case 'ode45'
                    if use_output_times
                        [t, Y] = ode45(f_rhs, tspan(:), v0, options);
                    else
                        [t, Y] = ode45(f_rhs, tspan, v0, options);
                    end
                case 'ode23'
                    if use_output_times
                        [t, Y] = ode23(f_rhs, tspan(:), v0, options);
                    else
                        [t, Y] = ode23(f_rhs, tspan, v0, options);
                    end
                otherwise
                    error('ODESolver: Unknown solver type: %s', obj.solver_type);
            end

            % Ensure Y is a normal matrix
            if isvector(Y) && size(Y,1) == 1
                Y = Y(:)';
            end
        end

        function setDefaultOptions(obj)
            % Base options only; AbsTol depends on state length so it is set in solve().
            obj.options = odeset('RelTol', 1.0e-12, 'Refine', 0);
        end

        function setOptions(obj, varargin)
            obj.options = odeset(obj.options, varargin{:});
        end

        function setSolver(obj, solver_type)
            valid_solvers = {'ode15s', 'ode23s', 'ode45', 'ode23'};
            if ~ismember(solver_type, valid_solvers)
                error('ODESolver: Invalid solver type. Choose from: %s', strjoin(valid_solvers, ', '));
            end
            obj.solver_type = solver_type;
        end

        function displayInfo(obj)
            fprintf('ODE Solver Configuration:\n');
            fprintf('  Solver type: %s\n', obj.solver_type);
        end
    end

    methods (Static)
        function [t, Y] = solveWithDefaults(rhs, tspan, v0)
            solver = ODESolver();
            [t, Y] = solver.solve(rhs, tspan, v0);
        end
    end
end