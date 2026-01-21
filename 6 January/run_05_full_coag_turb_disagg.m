% run_05_full_coag_turb_disagg.m
% Full physics: coag ON + turbulence ON + disagg ON (NEW disagg as a RATE)
% + (optional) explicit PP source if config supports it
%
% Key change vs your old version:
%   - Smooth epsilon pulse (solver-friendly, same experiment intent)
%   - Adds ode_options Stats + odeplot so you can SEE it's progressing

clear; close all; clc;

% ==========================================================
% 1) CONFIG
% ==========================================================
cfg = SimulationConfig();
cfg.use_column  = true;
cfg.z_max       = 65;
cfg.dz          = 5;
cfg.n_sections  = 24;

% -------------------------
% (Optional) Primary production source
% Only set if the property exists (prevents "Unrecognized property" errors)
% -------------------------
if isprop(cfg,'enable_pp'); cfg.enable_pp = true; end
if isprop(cfg,'pp_rate');   cfg.pp_rate   = 0.5; end   % example units: state per day
if isprop(cfg,'pp_bin');    cfg.pp_bin    = 1; end
if isprop(cfg,'pp_layer');  cfg.pp_layer  = 1; end

% -------------------------
% Physics switches
% -------------------------
if isprop(cfg,'enable_coag'); cfg.enable_coag = true; end

if isprop(cfg,'enable_disagg');     cfg.enable_disagg = true; end
if isprop(cfg,'disagg_apply_in');   cfg.disagg_apply_in = 'rhs'; end

% disaggregation RATE [d^-1]
if isprop(cfg,'disagg_rate'); cfg.disagg_rate = 1.0; end

% growth mode baseline
if isprop(cfg,'growth_mode'); cfg.growth_mode = 'shift'; end

% Optional RHS units debug
if isprop(cfg,'debug_rhs_units'); cfg.debug_rhs_units = false; end

% ==========================================================
% 2) TURBULENCE FORCING (smooth pulse)
% ==========================================================
eps0 = 1e-8;   % baseline
eps1 = 1e-6;   % pulse
t1   = 18;     % pulse start [d]
t2   = 20;     % pulse end [d]

% Smoothness timescale (days): 0.05–0.2 is typical
tau = 0.10;

if isprop(cfg,'eps_fun')
    cfg.eps_fun = @(t,z) ( eps0 + (eps1-eps0) .* ...
        (0.5*(1+tanh((t-t1)/tau))) .* (0.5*(1+tanh((t2-t)/tau))) );
end

% reference epsilon for disagg intensity MUST be baseline
if isprop(cfg,'eps_ref'); cfg.eps_ref = eps0; end

% ==========================================================
% 3) ODE OPTIONS (visibility + mild robustness)
% ==========================================================
% This does NOT change the physics, only helps you see progress and diagnose stiffness.
try
    cfg.ode_options = odeset( ...
        'Stats','on', ...          % prints solver stats at the end
        'OutputFcn', @odeplot ...  % live plot: shows it is moving
    );
catch
    % if your framework doesn't use cfg.ode_options, ignore
end

% ==========================================================
% 4) RUN SIMULATION
% ==========================================================
S   = CoagulationSimulation(cfg);
out = S.run();

% ==========================================================
% 5) DIAGNOSTICS
% ==========================================================
outdir = 'diag_05_full_coag_turb_disagg';
DiagnosticsSuite.runPart1(S, outdir);
DiagnosticsSuite.runPart2(S, outdir);

disp('RUN 05 DONE: full (coag + turb + disagg RATE)');

% ==========================================================
% 6) BUDGET CHECKS
% ==========================================================
debug_budget_one_time(S, 18);
debug_budget_one_time(S, 19);
debug_budget_one_time(S, 20);
debug_budget_one_time(S, 21);
debug_budget_one_time(S, 22);
debug_budget_one_time(S, 23);
debug_budget_one_time(S, 25);