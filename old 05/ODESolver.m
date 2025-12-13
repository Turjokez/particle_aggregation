classdef ODESolver < handle
    % ODESOLVER - Clean, robust adapter around MATLAB ODE solvers.
    % This version adds solve_odefun() for direct ODE function handles,
    % needed by the flux_test script (Adrian's mass/flux tests).

    properties
        solver_type = 'ode15s';
        options;
    end

    methods
        function obj = ODESolver(varargin)
            if nargin > 0
                obj.solver_type = varargin{1};
            end
            obj.setDefaultOptions();
        end

        % =====================================================================
        % SOLVE using CoagulationRHS object
        % =====================================================================
        function [t, Y] = solve(obj, rhs, tspan, v0, custom_options)
            if ~isa(rhs, 'CoagulationRHS')
                error('First argument to solve() must be a CoagulationRHS object');
            end
            rhs.validate();

            v0 = v0(:);
            n = numel(v0);

            % Merge options
            if nargin > 4 && ~isempty(custom_options)
                options = odeset(obj.options, custom_options);
            else
                options = obj.options;
            end

            % Non-negative constraint
            options = odeset(options, 'NonNegative', 1:n);

            % Jacobian (if implemented)
            try
                jacfun = @(tt,vv) rhs.jacobian(tt, vv(:));
                options = odeset(options, 'Jacobian', jacfun);
            catch
                % ignore if unavailable
            end

            % limit maximum internal time step
            dt = max(tspan(end) - tspan(1), eps);
            options = odeset(options, 'MaxStep', max(dt/4, 1e-6));

            % wrapper enforcing column shape + finite numbers
            f = @(tt,vv) local_rhs_wrapper(rhs, tt, vv);

            % Dispatch solver
            try
                [t, Y] = dispatch_solver(obj.solver_type, f, tspan, v0, options);
            catch ME
                % attempt interval split
                if dt <= 1e-6
                    rethrow(ME);
                end
                tm = 0.5*(tspan(1) + tspan(2));

                [t1, Y1] = dispatch_solver(obj.solver_type, f, [tspan(1) tm], v0, options);
                [t2, Y2] = dispatch_solver(obj.solver_type, f, [tm tspan(2)], Y1(end,:).', options);

                t = [t1; t2(2:end)];
                Y = [Y1; Y2(2:end,:)];
            end

            if size(Y,1)==1
                Y = Y.';
            end
        end

        % =====================================================================
        % NEW: Solve direct function handle (needed by flux_test)
        % =====================================================================
        function [t, Y] = solve_odefun(obj, fhandle, tspan, v0, custom_options)
            % fhandle(tt,yy) returns dv/dt
            v0 = v0(:);
            n  = numel(v0);

            % Merge options
            if nargin > 4 && ~isempty(custom_options)
                options = odeset(obj.options, custom_options);
            else
                options = obj.options;
            end

            % nonnegative constraint
            options = odeset(options, 'NonNegative', 1:n);

            % enforce MaxStep
            dt = max(tspan(end) - tspan(1), eps);
            options = odeset(options, 'MaxStep', max(dt/4, 1e-6));

            % wrap function for safety
            f = @(tt,vv) local_handle_wrapper(fhandle, tt, vv);

            % run solver
            [t, Y] = dispatch_solver(obj.solver_type, f, tspan, v0, options);

            if size(Y,1)==1
                Y = Y.';
            end
        end

        % =====================================================================
        function setDefaultOptions(obj)
            obj.options = odeset('RelTol', 1e-6, ...
                                 'AbsTol', 1e-12, ...
                                 'Refine', 0);
        end

        function setOptions(obj, varargin)
            obj.options = odeset(obj.options, varargin{:});
        end

        function setSolver(obj, solver_type)
            valid = {'ode15s','ode23s','ode45','ode23'};
            if ~ismember(solver_type, valid)
                error('Invalid solver type.');
            end
            obj.solver_type = solver_type;
        end

        function displayInfo(obj)
            fprintf('ODESolver type: %s\n', obj.solver_type);
            disp(obj.options);
        end
    end

    methods (Static)
        function [t, Y] = solveWithDefaults(rhs, tspan, v0)
            solver = ODESolver();
            [t, Y] = solver.solve(rhs, tspan, v0);
        end
    end
end

% ================================================================================
% LOCAL HELPERS
% ================================================================================

function dvdt = local_rhs_wrapper(rhs, tt, vv)
    dvdt = rhs.evaluate(tt, vv(:));
    dvdt = dvdt(:);
    dvdt(~isfinite(dvdt)) = 0;
end

function dvdt = local_handle_wrapper(fhandle, tt, vv)
    dvdt = fhandle(tt, vv(:));
    dvdt = dvdt(:);
    dvdt(~isfinite(dvdt)) = 0;
end

function [t, Y] = dispatch_solver(solver_type, f, tspan, v0, options)
    switch solver_type
        case 'ode15s'
            [t, Y] = ode15s(f, tspan, v0, options);
        case 'ode23s'
            [t, Y] = ode23s(f, tspan, v0, options);
        case 'ode45'
            [t, Y] = ode45(f, tspan, v0, options);
        case 'ode23'
            [t, Y] = ode23(f, tspan, v0, options);
        otherwise
            error('Unknown solver type %s', solver_type);
    end
end