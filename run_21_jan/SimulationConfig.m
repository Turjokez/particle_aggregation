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

        growth = 0.15;          % [d^{-1}]
        growth_mode = 'shift';  % 'shift' (legacy) or 'pp'

        gro_sec = 4;            % legacy shift start section
        num_1 = 1e3;            % legacy

        sinking_form = 'flux';  % 'flux' (recommended) or 'loop' (legacy)
        debug_sinking = false;

        % Disaggregation parameters (legacy kernel params)
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
        ode_options = [];
        clip_negative = true;
        use_full_jacobian = false;

        use_nonnegative = true;
        disagg_mode = 'legacy';

        enable_sinking = true;

        % ==========================================================
        % NEW-2025-12-13: coagulation on/off switches
        % ==========================================================
        enable_coag  = true;
        coag_scale   = 1;
        beta_scale   = 1;
        kernel_scale = 1;

        % Optional epsilon constant (RHS checks for this)
        epsilon_const = [];

        % ==========================================================
        % NEW-2025-12-20: turbulence-driven Disaggregation class controls
        % ==========================================================
        enable_disagg = false;
        disagg_apply_in = 'rhs';
        disagg_operator = 'applyWithScaling';

        eps_ref = 1e-6;
        disagg_n_exp = 0.45;

        disagg_kmax_a = 0.60;
        disagg_beta = 0.35;
        disagg_frac_to_edge = 1/3;
        disagg_redistribute_p = 0;
        enforce_coag_bv_closure = false;   % force coag to conserve biovolume by correcting last bin
        debug_coag_leak         = false;   % print leak info at runtime (optional)

        % Optional: time-varying epsilon forcing (for pulses)
        epsilon_time = [];
        epsilon_series = [];

        % ==========================================================
        % NEW-2025-12-21: depth-dependent turbulence forcing (eps(t,z))
        % ==========================================================
        eps_fun = [];
        epsilon_series_2d = [];
        epsilon_z = [];

        % ==========================================================
        % NEW-2025-12-21: choose ONE epsilon interface
        % ==========================================================
        epsilon_interface = 'auto';

        % ==========================================================
        % NEW-2025-12-21: RHS units debug printing (optional)
        % ==========================================================
        debug_rhs_units = false;

        % ==========================================================
        % NEW-2026-01-11: Disaggregation rate (per day) for RHS mode
        % ==========================================================
        disagg_rate = 1.0;

        % ==========================================================
        % NEW-2026-01-14: beta orientation switch (Step-2)
        % ==========================================================
        beta_row_form = 'legacy';
        beta_fix_b1_b3 = false;

        % ==========================================================
        % NEW-2026-01-15: explicit PP SOURCE term
        % ==========================================================
        enable_pp = false;
        pp_rate   = 0.0;
        pp_bin    = 1;
        pp_layer  = 1;

        % ==========================================================
        % NEW-2026-01-20: consistent export/inventory bin-weight mode
        % ==========================================================
        export_weight = "vbin";   % "ones" or "vbin"  <-- UPDATED DEFAULT
    end

    methods
        function obj = SimulationConfig(varargin)
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

                r = obj.z_max / obj.dz;
                if abs(r - round(r)) > 1e-12
                    warning('z_max/dz is not an integer (%.6g). getNumLayers() will floor(). Consider making z_max a multiple of dz for closure tests.', r);
                end
            end

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

            if ~isempty(obj.ode_options)
                assert(isstruct(obj.ode_options), 'ode_options must be an odeset struct or [].');
            end

            if isprop(obj,'debug_rhs_units') && ~isempty(obj.debug_rhs_units)
                assert(islogical(obj.debug_rhs_units) || isnumeric(obj.debug_rhs_units), ...
                    'debug_rhs_units must be true/false.');
            end

            % Disaggregation config checks
            if obj.enable_disagg
                validApply = any(strcmpi(obj.disagg_apply_in, {'simulation','rhs'}));
                assert(validApply, 'disagg_apply_in must be: simulation or rhs');

                validOp = any(strcmpi(obj.disagg_operator, {'applyWithScaling','apply'}));
                assert(validOp, 'disagg_operator must be: applyWithScaling or apply');

                assert(isfinite(obj.eps_ref) && obj.eps_ref > 0, 'eps_ref must be finite and > 0');
                assert(isfinite(obj.disagg_n_exp) && obj.disagg_n_exp >= 0, 'disagg_n_exp must be finite and >= 0');

                if isprop(obj,'disagg_rate') && ~isempty(obj.disagg_rate)
                    assert(isfinite(obj.disagg_rate) && obj.disagg_rate >= 0, 'disagg_rate must be finite and >= 0');
                end

                if ~isempty(obj.epsilon_time) || ~isempty(obj.epsilon_series)
                    assert(~isempty(obj.epsilon_time) && ~isempty(obj.epsilon_series), ...
                        'If using epsilon_time/epsilon_series, both must be provided.');
                end
            end

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
            % UPDATED: accept BOTH Nt-by-Nz and Nz-by-Nt for epsilon_series
            % ==========================================================
            if ~isempty(obj.epsilon_time) && ~isempty(obj.epsilon_series)

                % epsilon_time must be a vector (not 2D)
                assert(isnumeric(obj.epsilon_time) && isvector(obj.epsilon_time), ...
                    'epsilon_time must be a numeric VECTOR of times in days.');

                Nt = numel(obj.epsilon_time);

                if isvector(obj.epsilon_series)
                    assert(numel(obj.epsilon_series) == Nt, ...
                        'epsilon_series vector length must match epsilon_time length.');
                else
                    % matrix case
                    assert(isnumeric(obj.epsilon_series) && ismatrix(obj.epsilon_series), ...
                        'epsilon_series must be numeric (vector or matrix).');

                    if obj.use_column
                        Nz = obj.getNumLayers();

                        % Case 1: already Nt x Nz
                        if size(obj.epsilon_series,1) == Nt && size(obj.epsilon_series,2) == Nz
                            % ok
                        % Case 2: provided Nz x Nt -> transpose to Nt x Nz
                        elseif size(obj.epsilon_series,1) == Nz && size(obj.epsilon_series,2) == Nt
                            warning('epsilon_series provided as Nz-by-Nt. Transposing to Nt-by-Nz for internal use.');
                            obj.epsilon_series = obj.epsilon_series.'; % now Nt x Nz
                        else
                            error('epsilon_series matrix must be Nt-by-Nz OR Nz-by-Nt. Got %dx%d, expected Nt=%d, Nz=%d.', ...
                                size(obj.epsilon_series,1), size(obj.epsilon_series,2), Nt, Nz);
                        end
                    else
                        % 0-D: allow Nt x 1 or 1 x Nt or Nt x Nsomething? keep strict
                        if size(obj.epsilon_series,1) == Nt
                            % ok
                        elseif size(obj.epsilon_series,2) == Nt && size(obj.epsilon_series,1) == 1
                            obj.epsilon_series = obj.epsilon_series(:); % 1xNt -> Nt x 1
                        else
                            error('0-D epsilon_series must have Nt rows (or be a vector). Got %dx%d, Nt=%d.', ...
                                size(obj.epsilon_series,1), size(obj.epsilon_series,2), Nt);
                        end
                    end
                end
            end

            % ONE epsilon interface rule (warn on conflicts)
            has_fun    = ~isempty(obj.eps_fun);
            has_series = ~isempty(obj.epsilon_time) && ~isempty(obj.epsilon_series);
            has_2d     = ~isempty(obj.epsilon_series_2d) && ~isempty(obj.epsilon_z);
            has_const  = ~isempty(obj.epsilon_const);

            if strcmpi(string(obj.epsilon_interface),'auto')
                if (double(has_fun)+double(has_series)+double(has_2d)+double(has_const)) > 1
                    warning(['Multiple epsilon sources are set (eps_fun / epsilon_series / epsilon_series_2d / epsilon_const). ', ...
                             'RHS will use priority order. To avoid confusion, set cfg.epsilon_interface to one of: ', ...
                             '''fun'',''series'',''series_2d'',''const''.']);
                end
            else
                ei = lower(string(obj.epsilon_interface));
                if ei == "fun"
                    assert(has_fun, 'epsilon_interface="fun" but eps_fun is empty.');
                elseif ei == "series"
                    assert(has_series, 'epsilon_interface="series" but epsilon_time/epsilon_series are empty.');
                elseif ei == "series_2d"
                    assert(has_2d, 'epsilon_interface="series_2d" but epsilon_series_2d/epsilon_z are empty.');
                elseif ei == "const"
                    assert(has_const, 'epsilon_interface="const" but epsilon_const is empty.');
                else
                    error('epsilon_interface must be: auto|series|series_2d|fun|const');
                end
            end

            if isprop(obj,'beta_row_form') && ~isempty(obj.beta_row_form)
                bf = lower(string(obj.beta_row_form));
                validBF = any(bf == ["legacy","row_nrm","row_nr_m","row"]);
                assert(validBF, 'beta_row_form must be: legacy or row_nrm');
            end

            if isprop(obj,'beta_fix_b1_b3') && ~isempty(obj.beta_fix_b1_b3)
                assert(islogical(obj.beta_fix_b1_b3) || isnumeric(obj.beta_fix_b1_b3), ...
                    'beta_fix_b1_b3 must be true/false.');
            end

            if isprop(obj,'enable_pp') && ~isempty(obj.enable_pp) && logical(obj.enable_pp)
                assert(isfinite(obj.pp_rate) && obj.pp_rate >= 0, 'pp_rate must be finite and >= 0');
                assert(isfinite(obj.pp_bin) && obj.pp_bin >= 1 && obj.pp_bin <= obj.n_sections, ...
                    'pp_bin must be between 1 and n_sections');

                if obj.use_column
                    Nz = obj.getNumLayers();
                    assert(isfinite(obj.pp_layer) && obj.pp_layer >= 1 && obj.pp_layer <= Nz, ...
                        'pp_layer must be between 1 and Nz for column mode');
                end
            end

            % NEW: validate export_weight
            if isprop(obj,'export_weight') && ~isempty(obj.export_weight)
                ew = lower(string(obj.export_weight));
                assert(any(ew == ["ones","vbin"]), 'export_weight must be "ones" or "vbin".');
            end
        end

        function Nz = getNumLayers(obj)
            if obj.use_column
                Nz = max(1, floor(obj.z_max / obj.dz));
            else
                Nz = 1;
            end
        end

        function z = getZ(obj)
            Nz = obj.getNumLayers();
            if Nz == 1
                z = 0;
            else
                z = ((1:Nz) - 0.5) * obj.dz;
            end
        end
    end
end
