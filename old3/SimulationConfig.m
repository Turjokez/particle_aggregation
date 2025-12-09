classdef SimulationConfig < handle
    %SIMULATIONCONFIG Configuration for 1D Column Simulation
    
    properties
        % ===========================
        % Physical parameters
        % ===========================
        rho_fl      = 1.0275;      % Fluid density [g cm^{-3}]
        kvisc       = 0.01;        % Kinematic viscosity [cm^2 s^{-1}]
        g           = 980;         % Gravity [cm s^{-2}]
        day_to_sec  = 8.64e04;     % Seconds in a day
        k           = 1.3e-16;     % Boltzmann constant
        r_to_rg     = 1.36;        % Fractal scaling
        temp        = 20 + 273;    % Temp [K]
        salt        = 35;          % Salinity [psu]
        
        % ===========================
        % Section parameters
        % ===========================
        n_sections  = 20;          % Number of size bins
        kernel      = 'KernelBrown'; 
        d0          = 20e-4;       % Min diameter [cm]
        fr_dim      = 2.33;        
        num_1       = 10^3;        % Init concentration
        
        % ===========================
        % 1-D COLUMN GRID PARAMETERS
        % ===========================
        use_column (1,1) logical = true;   % Enable 1D mode
        z_max      = 200;                  % Max depth [m]
        dz         = 10;                   % Layer thickness [m]
        
        % ===========================
        % Biological & Physics parameters
        % ===========================
        alpha       = 1.0;         % Stickiness
        gamma       = 0.1;         % Shear rate [s^{-1}]
        growth      = 0.15;        % Specific growth rate [d^{-1}]
        gro_sec     = 4;           % Growth start section
        
        % NPP (Source)
        use_NPP         = false;
        NPP_rate        = 1e-4;
        NPP_profile     = 'constant';
        NPP_section     = 1;
        NPP_t_step      = 0;
        NPP_rate_after  = 0;
        NPP_amp         = 0.5;
        
        % ===========================
        % Legacy Linear Breakup (Restored)
        % ===========================
        c3 = 0.1;   % Fragmentation rate const
        c4 = 1.4;   % Fragment distribution slope
        
        % ===========================
        % Turbulence & Disaggregation (New)
        % ===========================
        epsilon_profile char = 'observed'; 
        epsilon          = 1e-7;
        epsilon_time     = [];     
        epsilon_series   = [];
        epsilon_ref      = 1e-6; 
        epsilon_mean     = 1e-6;
        epsilon_amp      = 0;
        epsilon_period   = 30;
        epsilon_phase    = 0;
        
        % Scaling of alpha with epsilon
        alpha_base       = 1.0;
        p_alpha          = 0.0; 
        alpha_clip_min   = 1e-3;
        alpha_clip_max   = 50;

        % Nonlinear Breakup Settings
        disagg_use_nonlinear (1,1) logical = true;
        disagg_kmax_a           = 0.70; 
        disagg_beta             = 0.35;
        disagg_frac_to_edge     = 0.40;
        disagg_redistribute_p   = 0.50;
        
        % ===========================
        % Time Integration
        % ===========================
        t_init  = 0.0;
        t_final = 30.0;
        delta_t = 0.5;
        tracer  = false;
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
            if mod(obj.z_max, obj.dz) > 1e-5
                warning('z_max is not a multiple of dz. Grid may be uneven.');
            end
        end
    end
end