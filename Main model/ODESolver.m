classdef ODESolver < handle
    %ODESOLVER Adapter for MATLAB ODE solvers

    properties
        solver_type = 'ode15s';  % Default solver
        options;                 % ODE solver options
    end

    methods
        function obj = ODESolver(varargin)
            %ODESOLVER Constructor
            if nargin > 0
                obj.solver_type = varargin{1};
            end

            % Set default options
            obj.setDefaultOptions();
        end

        function [t, Y] = solve(obj, rhs, tspan, v0, custom_options)
            %SOLVE Solve ODE system
            % rhs = CoagulationRHS object
            % tspan = time span vector
            % v0 = initial conditions vector
            % custom_options = optional custom ODE options

            % Validate inputs
            if ~isa(rhs, 'CoagulationRHS')
                error('First argument must be a CoagulationRHS object');
            end

            rhs.validate();

            % Get number of sections from initial conditions
            n_sections = length(v0);

            % Set up adaptive absolute tolerance based on actual number of sections
            abs_tol = 1.0e-18;
            calcomp = 1:n_sections;
            at = (abs_tol * 1.5 .^ (-(calcomp-1)));

            % Create options with correct AbsTol and allow negative values
            options = odeset('RelTol', 3.0e-14, 'Refine', 0, 'AbsTol', at, 'NonNegative', []);

            % Merge custom options if provided
            if nargin > 4 && ~isempty(custom_options)
                options = odeset(options, custom_options);
            end

            % Set up Jacobian function handle
            jacobian_fun = @(t, v) rhs.jacobian(t, v);
            options = odeset(options, 'Jacobian', jacobian_fun);

            % Solve ODE
            switch obj.solver_type
                case 'ode15s'
                    [t, Y] = ode15s(@(t, v) rhs.evaluate(t, v), tspan, v0, options);
                case 'ode23s'
                    [t, Y] = ode23s(@(t, v) rhs.evaluate(t, v), tspan, v0, options);
                case 'ode45'
                    [t, Y] = ode45(@(t, v) rhs.evaluate(t, v), tspan, v0, options);
                case 'ode23'
                    [t, Y] = ode23(@(t, v) rhs.evaluate(t, v), tspan, v0, options);
                otherwise
                    error('Unknown solver type: %s', obj.solver_type);
            end

            % Ensure solution is a matrix (ode15s can return column vectors)
            if size(Y, 1) == 1
                Y = Y';
            end
        end

        function setDefaultOptions(obj)
            %SETDEFAULTOPTIONS Set default ODE solver options
            n_sections = 20;  % Default, will be overridden
            abs_tol = 1.0e-18;
            rel_tol = 3.0e-14;

            % Adaptive absolute tolerance
            calcomp = 1:n_sections;
            at = (abs_tol * 1.5 .^ (-(calcomp-1)));

            obj.options = odeset('RelTol', rel_tol, 'Refine', 0, 'AbsTol', at, 'NonNegative', []);
        end

        function setOptions(obj, varargin)
            %SETOPTIONS Set ODE solver options using odeset syntax
            obj.options = odeset(obj.options, varargin{:});
        end

        function setSolver(obj, solver_type)
            %SETSOLVER Change solver type
            valid_solvers = {'ode15s', 'ode23s', 'ode45', 'ode23'};
            if ~ismember(solver_type, valid_solvers)
                error('Invalid solver type. Choose from: %s', strjoin(valid_solvers, ', '));
            end
            obj.solver_type = solver_type;
        end

        function displayInfo(obj)
            %DISPLAYINFO Display solver configuration
            fprintf('ODE Solver Configuration:\n');
            fprintf('  Solver type: %s\n', obj.solver_type);
            fprintf('  Options:\n');

            % Display key options
            if isfield(obj.options, 'RelTol')
                fprintf('    Relative tolerance: %.2e\n', obj.options.RelTol);
            end
            if isfield(obj.options, 'AbsTol')
                if isscalar(obj.options.AbsTol)
                    fprintf('    Absolute tolerance: %.2e\n', obj.options.AbsTol);
                else
                    fprintf('    Absolute tolerance: vector [%.2e, %.2e, ...]\n', ...
                        obj.options.AbsTol(1), obj.options.AbsTol(end));
                end
            end
            if isfield(obj.options, 'Refine')
                fprintf('    Refine: %d\n', obj.options.Refine);
            end
        end
    end

    methods (Static)
        function [t, Y] = solveWithDefaults(rhs, tspan, v0)
            %SOLVEWITHDEFAULTS Solve ODE with default solver and options
            % Convenience method for quick solving

            solver = ODESolver();
            [t, Y] = solver.solve(rhs, tspan, v0);
        end
    end
end
