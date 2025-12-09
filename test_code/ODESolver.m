classdef ODESolver < handle
    %ODESOLVER Adapter for MATLAB ODE solvers (robust, column-safe)

    properties
        solver_type = 'ode15s';  % Default solver
        options;                 % Default ODE solver options
    end

    methods
        function obj = ODESolver(varargin)
            % ODESOLVER Constructor
            if nargin > 0
                obj.solver_type = varargin{1};
            end
            obj.setDefaultOptions();
        end

        function [t, Y] = solve(obj, rhs, tspan, v0, custom_options)
            % SOLVE Solve ODE system
            %   rhs           = CoagulationRHS object
            %   tspan         = [t0 tf]
            %   v0            = initial state (row/col ok; will be forced col)
            %   custom_options= optional odeset(options)

            if ~isa(rhs, 'CoagulationRHS')
                error('First argument must be a CoagulationRHS object');
            end
            rhs.validate();

            % --- ensure column state and build per-interval options ---
            v0 = v0(:);
            n_sections = numel(v0);

            % merge options: class defaults -> caller overrides
            if nargin > 4 && ~isempty(custom_options)
                options = odeset(obj.options, custom_options);
            else
                options = obj.options;
            end

            % NonNegative on all species
            options = odeset(options, 'NonNegative', 1:n_sections);

            % attach Jacobian (shape-safe)
            try
                jacfun  = @(tt,vv) rhs.jacobian(tt, vv(:));
                options = odeset(options, 'Jacobian', jacfun);
            catch
                % ok to skip if rhs has no jacobian
            end

            % cap internal step to help around sharp forcing
            dt = max(tspan(end) - tspan(1), eps);
            options = odeset(options, 'MaxStep', max(dt/4, 1e-6));  % time units = days

            % RHS wrapper guarantees column output & finite numbers
            f = @(tt,vv) local_rhs(rhs, tt, vv);

            % --- try chosen solver; on failure, split interval and recurse ---
            try
                [t, Y] = dispatch_solver(obj.solver_type, f, tspan, v0, options);
            catch ME
                if dt <= 1e-6
                    rethrow(ME);
                end
                tm = 0.5*(tspan(1) + tspan(2));

                [t1, Y1] = dispatch_solver(obj.solver_type, f, [tspan(1) tm], v0, options);
                v_mid     = Y1(end,:).';  % column

                [t2, Y2] = dispatch_solver(obj.solver_type, f, [tm tspan(2)], v_mid, options);

                % concatenate without duplicating the midpoint
                t = [t1; t2(2:end)];
                Y = [Y1; Y2(2:end,:)];
            end

            % ensure Y is [nt x n] even if solver returns a single row
            if size(Y,1) == 1
                Y = Y.';
            end
        end

        function setDefaultOptions(obj)
            % SETDEFAULTOPTIONS: sane, stiff-friendly defaults
            n_sections = 20;        % placeholder; per-call NonNegative is set in solve
            rel_tol    = 1e-6;      % much looser than 1e-14 => avoids stiffness dead-ends
            abs_tol    = 1e-12;     % scalar; robust baseline

            % Keep vector AbsTol style if you prefer — leave scalar by default:
            % calcomp = 1:n_sections;
            % at = (abs_tol * 1.5 .^ (-(calcomp-1)));

            obj.options = odeset('RelTol', rel_tol, ...
                                 'AbsTol', abs_tol, ...
                                 'Refine', 0);   % keep dense-output off
        end

        function setOptions(obj, varargin)
            % SETOPTIONS Set ODE solver options using odeset syntax
            obj.options = odeset(obj.options, varargin{:});
        end

        function setSolver(obj, solver_type)
            % SETSOLVER Change solver type
            valid_solvers = {'ode15s','ode23s','ode45','ode23'};
            if ~ismember(solver_type, valid_solvers)
                error('Invalid solver type. Choose from: %s', strjoin(valid_solvers, ', '));
            end
            obj.solver_type = solver_type;
        end

        function displayInfo(obj)
            % DISPLAYINFO Display solver configuration
            fprintf('ODE Solver Configuration:\n');
            fprintf('  Solver type: %s\n', obj.solver_type);
            s = obj.options;
            if ~isempty(s)
                if isfield(s,'RelTol') && ~isempty(s.RelTol)
                    fprintf('  RelTol: %.2e\n', s.RelTol);
                end
                if isfield(s,'AbsTol') && ~isempty(s.AbsTol)
                    if isscalar(s.AbsTol)
                        fprintf('  AbsTol: %.2e\n', s.AbsTol);
                    else
                        fprintf('  AbsTol: vector [%g ... %g]\n', s.AbsTol(1), s.AbsTol(end));
                    end
                end
                if isfield(s,'MaxStep') && ~isempty(s.MaxStep)
                    fprintf('  MaxStep: %g\n', s.MaxStep);
                end
                if isfield(s,'NonNegative') && ~isempty(s.NonNegative)
                    fprintf('  NonNegative set on %d variables\n', numel(s.NonNegative));
                end
            end
        end
    end

    methods (Static)
        function [t, Y] = solveWithDefaults(rhs, tspan, v0)
            % SOLVEWITHDEFAULTS Convenience solve with defaults
            solver = ODESolver();
            [t, Y] = solver.solve(rhs, tspan, v0);
        end
    end
end

% ===== local helper functions (not methods) ===============================

function dvdt = local_rhs(rhs, tt, vv)
    % Guarantee column in/out and finite values for the integrator
    dvdt = rhs.evaluate(tt, vv(:));
    dvdt = dvdt(:);
    dvdt(~isfinite(dvdt)) = 0;
end

function [t, Y] = dispatch_solver(solver_type, f, tspan, v0, options)
    switch solver_type
        case 'ode15s'
            [t, Y] = ode15s(f, tspan, v0(:), options);
        case 'ode23s'
            [t, Y] = ode23s(f, tspan, v0(:), options);
        case 'ode45'
            [t, Y] = ode45( f, tspan, v0(:), options);
        case 'ode23'
            [t, Y] = ode23( f, tspan, v0(:), options);
        otherwise
            error('Unknown solver type: %s', solver_type);
    end
end