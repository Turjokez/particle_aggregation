classdef SimulationConfig < handle
    %SIMULATIONCONFIG Configuration class for coagulation simulation parameters
    
    properties
        % Physical parameters
        rho_fl = 1.0275;        % Fluid density [g cm^{-3}]
        kvisc = 0.01;           % Kinematic viscosity [cm^2 s^{-1}]
        g = 980;                % Accel. due to gravity [cm s^{-2}]
        day_to_sec = 8.64e04;   % Seconds in a day [s d^{-1}]
        k = 1.3e-16;           % Boltzmanns constant [erg K^{-1}]
        r_to_rg = 1.36;        % Interaction radius to radius of gyration
        
        % Section/coagulation related parameters
        n_sections = 20;        % Number of sections
        kernel = 'KernelBrown'; % Kernel type
        d0 = 20e-4;            % Diameter of unit particle [cm]
        fr_dim = 2.33;         % Particle fractal dimension
        n1 = 100;              % No. particles cm^{-3} in first section
        
        % Other input parameters
        temp = 20 + 273;       % Temperature [K]
        alpha = 1.0;           % Stickiness
        dz = 65;               % Layer thickness [m]
        gamma = 0.1;           % Average shear rate [s^{-1}]
        growth = 0.15;         % Specific growth rate in first section [d^{-1}]
        gro_sec = 4;           % Section at which growth in aggregates starts
        num_1 = 10^3;          % Number of particle cm^{-3} in first section
        
        % Disaggregation parameters
        c3 = 0.2;              % For curvilinear kernel
        c4 = 1.45;             % For curvilinear kernel
        
        % Parameters for solving equations
        t_init = 0.0;          % Initial time for integrations [d]
        t_final = 30.0;        % Final time for integrations [d]
        delta_t = 1.0;         % Time interval for output [d]
        
        % Code Runtime Options
        tracer = false;         % Integrate tracer as well [false=no, true=yes]
    end
    
    methods
        function obj = SimulationConfig(varargin)
            % Constructor - can accept parameter-value pairs
            if nargin > 0
                for i = 1:2:length(varargin)
                    if isprop(obj, varargin{i})
                        obj.(varargin{i}) = varargin{i+1};
                    end
                end
            end
        end
        
        function grid = derive(obj)
            %DERIVE Compute derived grid and constants
            % Returns: DerivedGrid object with precomputed values
            grid = DerivedGrid(obj);
        end
        
        function validate(obj)
            %VALIDATE Validate configuration parameters
            assert(obj.n_sections > 0, 'n_sections must be positive');
            assert(obj.t_final > obj.t_init, 't_final must be greater than t_init');
            assert(obj.delta_t > 0, 'delta_t must be positive');
            % Add more validation as needed
        end
    end
end
