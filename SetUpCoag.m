function [param, opt] = SetUpCoag()
%
% SetUpCoag obtains user options for the coagulation calculations and then
% does some housecleaning
%
% USAGE:
%   
%
% HISTORY:
%   23-04-09: First cut.
%
%
% Adrian Burd, University of Georgia, 2009
%

%% Physical parameters

param.rho_fl     = 1.0275;        % Fluid density [g cm^{-3}]
param.kvisc      = 0.01;          % Kinematic viscosity [cm^2 s^{-1}]
param.g          = 980;           % Accel. due to gravity [cm s^{-2}]
param.day_to_sec = 8.64e04;       % Seconds in a day [s d^{-1}]
param.k          = 1.3e-16;       % Boltzmanns constant [erg K^{-1}]
param.r_to_rg    = 1.36;          % Interaction radius to radius of gyration

%% Section/coagulation related parameters

param.n_sections = 20;            % Nmber of sections
param.kernel     = 'KernelBrown'; % Kernel type
param.d0         = 20e-4;         % Diameter of unit particle [cm] [20e-4]
param.fr_dim     = 2.33;          % Particle fractal dimension
param.n1         = 100;           % No. particles cm^{-3} in first section

%% Other input parameters - set up input at a later date

param.temp    = 20 + 273;         % Temperature [K]
param.alpha   = 1.0;              % Stickiness
param.dz      = 65;               % Layer thickness [m]
param.gamma   = 0.1;              % Average shear rate [s^{-1}]
param.growth  = 0.15;             % Specific growth rate in first section [d^{-1}]
param.gro_sec = 4;                % Section at which growth in aggregates starts
param.num_1   = 10^3;               % Number of particle cm^{-3} in first section [40]

%% Disaggregation parameters
%
% FOr curvilinear kernel, c3=0.20, c4=1.45. For rectliniear kernel c3=13.0,
% c4=1.31.

param.c3 = 0.2;
param.c4 = 1.45;

%% Paremeters for solving equations

param.t_init  = 0.0;              % Initial time for integrations [d]
param.t_final = 30.0;             % Final time for integrations [d]
param.delta_t = 1.0;              % Time interval for output [d]

%% Code Runtime Options

opt.tracer = 0;                   % Integrate tracer as well [0=no, 1=yes]

%% Derived parameters

param.dvisc = param.kvisc*param.rho_fl;   % Dynamic viscosity [g cm^{-1} s^{-1}]
param.del_rho = (4.5*2.48)*param.kvisc*param.rho_fl/param.g*(param.d0/2)^(-0.83);

param.conBr = 2.0/3.0*param.k*param.temp/param.dvisc;

param.a0 = param.d0/2;
param.v0 = (pi/6) * param.d0^3;

param.v_lower = param.v0*2.^( 0 : param.n_sections-1)';
param.v_upper = 2.0*param.v_lower;

param.av_vol  = 1.5*param.v_lower;
param.dcomb   = (param.v_lower*6/pi).^(1.0/3.0);
param.dwidth  = (2^(1.0/3.0)-1)*param.dcomb;

amfrac       = (4.0/3.0*pi)^(-1.0/param.fr_dim) * param.a0^(1.0-3.0/param.fr_dim);
param.amfrac = amfrac*sqrt(0.6);
param.bmfrac = 1./param.fr_dim;

param.setcon = (2.0/9.0)*param.del_rho/param.rho_fl*param.g/param.kvisc;


%% Tracer parameters
%
% Store parameters for the particular trace metal or tracer here
%


