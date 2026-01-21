classdef SimulationConfig < handle
    %SIMULATIONCONFIG Configuration class for coagulation simulation parameters

    properties
        % Physical parameters
        rho_fl = 1.0275;        % Fluid density [g cm^{-3}]
        kvisc = 0.01;           % Kinematic viscosity [cm^2 s^{-1}]
        g = 980;                % Accel. due to gravity [cm s^{-2}]
        day_to_sec = 8.64e04;   % Seconds in a day [s d^{-1}]
        k = 1.3e-16;            % Boltzmanns constant [erg K^{-1}]
        r_to_rg = 1.36;         % Interaction radius to radius of gyration

        % Section/coagulation related parameters
        n_sections = 20;        % Number of sections
        kernel = 'KernelBrown'; % Kernel type
        d0 = 20e-4;             % Diameter of unit particle [cm]
        fr_dim = 2.33;          % Particle fractal dimension
        n1 = 100;               % No. particles cm^{-3} in first section

        % Other input parameters
        temp = 20 + 273;        % Temperature [K]
        alpha = 1.0;            % Stickiness

        % ---- dz meaning ----
        % For 0-D slab: dz is the slab thickness [m]
        % For 1-D column: dz is the layer thickness [m]
        dz = 65;                % [m] default slab thickness

        % ---- Column options ----
        use_column = false;     % false=0-D slab, true=1-D column
        z_max      = 65;        % total depth [m]

        % ---- Initial condition vertical profile options ----
        init_profile = '';      % 'surface_only'|'uniform'|'exp'|'top_only' or '' for default
        init_z0 = 30;           % e-fold depth (m) if init_profile='exp'

        % ---- Linear processes ----
        gamma = 0.1;            % Average shear rate [s^{-1}]

        % growth is the "rate parameter" used by LinearProcessBuilder
        %  - in legacy shift-mode it is the legacy growth strength
        %  - in pp-mode it is mu (1/d)
        growth = 0.15;          % [d^{-1}]
        growth_mode = 'shift';  % 'shift' (legacy) or 'pp' (Adrian Test 1)

        gro_sec = 4;            % legacy shift start section
        num_1 = 1e3;            % legacy

        % Optional: choose sinking operator form
        sinking_form = 'flux';  % 'flux' (recommended) or 'loop' (legacy)

        % Optional: debug switch for sinking conservation check
        debug_sinking = false;

        % Disaggregation parameters
        c3 = 0.2;               % For curvilinear kernel
        c4 = 1.45;              % For curvilinear kernel

        % Time integration
        t_init = 0.0;           % [d]
        t_final = 30.0;         % [d]
        delta_t = 1.0;          % [d]

        % Code Runtime Options
        tracer = false;         % Integrate tracer as well [false=no, true=yes]

        % ==========================================================
        % NEW (safe defaults)
        % ==========================================================
        ode_options = [];          % NEW: allows cfg.ode_options without error (odeset struct or [])
        clip_negative = true;      % NEW: clamp negatives to 0 inside RHS
        use_full_jacobian = false; % NEW: keep Jacobian safe unless you explicitly enable

        % ----------------------------------------------------------
        % NEW: allow turning solver NonNegative on/off for budget tests
        % ----------------------------------------------------------
        use_nonnegative = true;    % NEW
        disagg_mode = 'legacy';    % or 'loop_equiv'

        % ==========================================================
        % NEW-2025-12-13: coagulation on/off switches (needed for tests 2A/2B/2C)
        % ==========================================================
        enable_coag  = true;   % master switch (false => skip kernel compute)
        coag_scale   = 1;      % 0 disables coag RHS contribution
        beta_scale   = 1;      % 0 disables betas
        kernel_scale = 1;      % 0 disables kernel compute path

        % Optional epsilon constant (RHS checks for this)
        epsilon_const = [];        % NEW: fallback epsilon (only stored unless you use it in kernels)

        % ==========================================================
        % NEW-2025-12-20: turbulence-driven Disaggregation class controls
        % (used by Disaggregation.applyWithScaling)
        % ==========================================================
        enable_disagg = false;              % NEW-2025-12-20: master switch for Disaggregation.m

        % NOTE (NEW-2025-12-20):
        %   Use 'rhs' to apply Disaggregation.applyWithScaling inside CoagulationRHS
        %   and zero out matrix disaggregation in CoagulationSimulation (prevents double counting).
        %   Use 'simulation' ONLY if you plan to apply breakup outside RHS (not currently used).
        disagg_apply_in = 'rhs';            % NEW-2025-12-20: 'rhs' (recommended) or 'simulation' (advanced)

        disagg_operator = 'applyWithScaling'; % NEW-2025-12-20: 'applyWithScaling' or 'apply' (legacy API)

        % epsilon control for disaggregation scaling
        eps_ref = 1e-6;                     % NEW-2025-12-20: reference ε [W kg^-1]
        disagg_n_exp = 0.45;                % NEW-2025-12-20: nonlinearity exponent for intensity

        % Optional tuning parameters used by Disaggregation.applyWithScaling
        disagg_kmax_a = 0.60;               % NEW-2025-12-20
        disagg_beta = 0.35;                 % NEW-2025-12-20
        disagg_frac_to_edge = 1/3;          % NEW-2025-12-20
        disagg_redistribute_p = 0;          % NEW-2025-12-20

        % Optional: time-varying epsilon forcing (for pulses)
        epsilon_time = [];                  % NEW-2025-12-20: vector of times [d]
        epsilon_series = [];                % NEW-2025-12-20: vector (0-D) OR Nt-by-Nz (column) OR 1-by-Nt

        % ==========================================================
        % NEW-2025-12-21: depth-dependent turbulence forcing (eps(t,z))
        % Priority recommendation:
        %   1) eps_fun(t_days, z_m)
        % Optional storage:
        %   epsilon_series_2d (Nt x Nz) + epsilon_z (Nz x 1)
        % ==========================================================
        eps_fun = [];              % NEW-2025-12-21: function handle: eps_fun(t, z) -> epsilon [W kg^-1]
        epsilon_series_2d = [];    % NEW-2025-12-21: Nt x Nz epsilon(t,z)
        epsilon_z = [];            % NEW-2025-12-21: Nz x 1 depth grid (m) for epsilon_series_2d

        % ==========================================================
        % NEW-2025-12-21: choose ONE epsilon interface
        % ==========================================================
        epsilon_interface = 'auto'; % NEW-2025-12-21: 'auto'|'series'|'series_2d'|'fun'|'const'

        % ==========================================================
        % NEW-2025-12-21: RHS units debug printing (optional)
        % Used by CoagulationRHS to print ||A*v|| and ||rhs|| ratio
        % ==========================================================
        debug_rhs_units
        = false;    % NEW-2025-12-21
    end

    methods
        function obj = SimulationConfig(varargin)
            % Constructor: parameter-value pairs
            if nargin > 0
                for i = 1:2:length(varargin)
                    if isprop(obj, varargin{i})
                        obj.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end

        function grid = derive(obj)
            grid = DerivedGrid(obj);
        end

        function validate(obj)
            assert(obj.n_sections > 0, 'n_sections must be positive');
            assert(obj.t_final > obj.t_init, 't_final must be greater than t_init');
            assert(obj.delta_t > 0, 'delta_t must be positive');

            if obj.use_column
                assert(obj.dz > 0, 'dz must be positive when use_column = true');
                assert(obj.z_max > 0, 'z_max must be positive when use_column = true');
                assert(obj.z_max >= obj.dz, 'z_max must be >= dz when use_column = true');

                % NEW-2025-12-12: warn if z_max not divisible by dz (can change closure/export)
                r = obj.z_max / obj.dz;
                if abs(r - round(r)) > 1e-12
                    warning('z_max/dz is not an integer (%.6g). getNumLayers() will floor(). Consider making z_max a multiple of dz for closure tests.', r);
                end
            end

            % ==========================================================
            % NEW-2025-12-12: allow 'top_only' for Adrian Test 2
            % ==========================================================
            if obj.use_column && ~isempty(obj.init_profile)
                valid = any(strcmpi(obj.init_profile, {'surface_only','uniform','exp','top_only'}));
                assert(valid, 'init_profile must be: surface_only, uniform, exp, or top_only');

                if strcmpi(obj.init_profile,'exp')
                    assert(obj.init_z0 > 0, 'init_z0 must be positive for exp profile');
                end
            end

            if ~isempty(obj.growth_mode)
                validGM = any(strcmpi(obj.growth_mode, {'shift','pp'}));
                assert(validGM, 'growth_mode must be: shift or pp');
            end

            if ~isempty(obj.sinking_form)
                validSF = any(strcmpi(obj.sinking_form, {'flux','loop'}));
                assert(validSF, 'sinking_form must be: flux or loop');
            end

            % ==========================================================
            % NEW: ode_options sanity check
            % ==========================================================
            if ~isempty(obj.ode_options)
                assert(isstruct(obj.ode_options), 'ode_options must be an odeset struct or [].');
            end

            % ==========================================================
            % NEW-2025-12-21: debug_rhs_units sanity check (safe)
            % ==========================================================
            if isprop(obj,'debug_rhs_units') && ~isempty(obj.debug_rhs_units)
                assert(islogical(obj.debug_rhs_units) || isnumeric(obj.debug_rhs_units), ...
                    'debug_rhs_units must be true/false.');
            end

            % ==========================================================
            % NEW-2025-12-20: Disaggregation config checks
            % ==========================================================
            if obj.enable_disagg
                validApply = any(strcmpi(obj.disagg_apply_in, {'simulation','rhs'}));
                assert(validApply, 'disagg_apply_in must be: simulation or rhs');

                validOp = any(strcmpi(obj.disagg_operator, {'applyWithScaling','apply'}));
                assert(validOp, 'disagg_operator must be: applyWithScaling or apply');

                assert(isfinite(obj.eps_ref) && obj.eps_ref > 0, 'eps_ref must be finite and > 0');
                assert(isfinite(obj.disagg_n_exp) && obj.disagg_n_exp >= 0, 'disagg_n_exp must be finite and >= 0');

                % If user provides epsilon time series, basic consistency check
                if ~isempty(obj.epsilon_time) || ~isempty(obj.epsilon_series)
                    assert(~isempty(obj.epsilon_time) && ~isempty(obj.epsilon_series), ...
                        'If using epsilon_time/epsilon_series, both must be provided.');
                end
            end

            % ==========================================================
            % NEW-2025-12-21: depth-dependent turbulence forcing checks
            % ==========================================================
            if ~isempty(obj.eps_fun)
                assert(isa(obj.eps_fun,'function_handle'), 'eps_fun must be a function handle: eps_fun(t,z)->epsilon.');
            end

            if ~isempty(obj.epsilon_series_2d) || ~isempty(obj.epsilon_z)
                assert(~isempty(obj.epsilon_series_2d) && ~isempty(obj.epsilon_z), ...
                    'If using epsilon_series_2d/epsilon_z, both must be provided.');
                assert(isnumeric(obj.epsilon_series_2d) && ismatrix(obj.epsilon_series_2d), ...
                    'epsilon_series_2d must be numeric Nt x Nz.');
                assert(isnumeric(obj.epsilon_z) && isvector(obj.epsilon_z), ...
                    'epsilon_z must be a numeric vector (Nz x 1).');
            end

            % ==========================================================
            % NEW-2025-12-21: strengthen epsilon_series shape expectations
            % epsilon_series is allowed to be:
            %   - Nt x 1  (vector, old)
            %   - 1 x Nt  (row vector)
            %   - Nt x Nz (column depth dependent)
            % ==========================================================
            if ~isempty(obj.epsilon_time) && ~isempty(obj.epsilon_series)
                Nt = numel(obj.epsilon_time);

                if isvector(obj.epsilon_series)
                    assert(numel(obj.epsilon_series) == Nt, ...
                        'epsilon_series vector length must match epsilon_time length.');
                else
                    % matrix case
                    assert(size(obj.epsilon_series,1) == Nt, ...
                        'epsilon_series matrix must be Nt-by-Nz (rows=time).');
                    if obj.use_column
                        Nz = obj.getNumLayers();
                        assert(size(obj.epsilon_series,2) == Nz, ...
                            'epsilon_series matrix must have Nz columns to match config grid.');
                    end
                end
            end

            % ==========================================================
            % NEW-2025-12-21: ONE epsilon interface rule (prevents silent conflicts)
            % ==========================================================
            has_fun    = ~isempty(obj.eps_fun);
            has_series = ~isempty(obj.epsilon_time) && ~isempty(obj.epsilon_series);
            has_2d     = ~isempty(obj.epsilon_series_2d) && ~isempty(obj.epsilon_z);
            has_const  = ~isempty(obj.epsilon_const);

            % OLD (kept): no rule existed, could provide multiple at once.
            % NEW: enforce that user does not accidentally provide multiple sources
            if strcmpi(string(obj.epsilon_interface),'auto')
                % allow multiple to exist, but we warn (RHS uses a priority order)
                if (double(has_fun)+double(has_series)+double(has_2d)+double(has_const)) > 1
                    warning(['Multiple epsilon sources are set (eps_fun / epsilon_series / epsilon_series_2d / epsilon_const). ', ...
                             'RHS will use priority order. To avoid confusion, set cfg.epsilon_interface to one of: ', ...
                             '''fun'',''series'',''series_2d'',''const''.']);
                end
            else
                % strict mode: only one is allowed depending on epsilon_interface
                ei = lower(string(obj.epsilon_interface));
                if ei == "fun"
                    assert(has_fun, 'epsilon_interface="fun" but eps_fun is empty.');
                elseif ei == "series"
                    assert(has_series, 'epsilon_interface="series" but epsilon_time/epsilon_series are empty.');
                elseif ei == "series_2d"
                    assert(has_2d, 'epsilon_interface="series_2d" but epsilon_series_2d/epsilon_z are empty.');
                elseif ei == "const"
                    assert(has_const, 'epsilon_interface="const" but epsilon_const is empty.');
                end
            end
        end

        function Nz = getNumLayers(obj)
            % For 0-D slab, returns 1.
            if obj.use_column
                % Keep floor() so grid is stable and predictable.
                % (If you want exact depth closure, make z_max a multiple of dz.)
                Nz = max(1, floor(obj.z_max / obj.dz));
            else
                Nz = 1;
            end
        end

        function z = getZ(obj)
            % Depth grid (m), cell centers
            Nz = obj.getNumLayers();
            if Nz == 1
                z = 0;
            else
                z = ((1:Nz) - 0.5) * obj.dz;
            end
        end
    end
end