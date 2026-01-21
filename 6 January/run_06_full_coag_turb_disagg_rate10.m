% run_06_full_coag_turb_disagg_rate10.m
% Goal: Full physics: coag ON + turbulence ON + disagg ON (NEW disagg as a RATE)

clear; close all;

cfg = SimulationConfig();
cfg.use_column  = true;
cfg.z_max       = 65;
cfg.dz          = 5;
cfg.n_sections  = 24;

% -------------------------
% Physics switches
% -------------------------
cfg.enable_coag = true;

if isprop(cfg,'enable_disagg');     cfg.enable_disagg = true; end
if isprop(cfg,'disagg_apply_in');   cfg.disagg_apply_in = 'rhs'; end

% ---------------------------------------------------------
% NEW-2026-01-11: disaggregation RATE (per day)
% ---------------------------------------------------------
if isprop(cfg,'disagg_rate')
    cfg.disagg_rate = 10.0;   % [d^-1]
end

% Keep baseline
cfg.growth_mode = 'shift';

% -------------------------
% Turbulence forcing: pulse (UPDATED TIME)
% -------------------------
eps0 = 1e-8;
eps1 = 1e-6;
t1   = 18;      % <<< UPDATED
t2   = 20;      % <<< UPDATED

if isprop(cfg,'eps_fun')
    cfg.eps_fun = @(t,z) (eps0 + (eps1-eps0) * double(t >= t1 & t <= t2));
end

% eps_ref should be your baseline (so eps_rel = eps/eps_ref)
if isprop(cfg,'eps_ref')
    cfg.eps_ref = eps0;
end

% Optional (turn OFF if you don't want prints)
if isprop(cfg,'debug_rhs_units')
    cfg.debug_rhs_units = false;
end

S = CoagulationSimulation(cfg);
out = S.run();

DiagnosticsSuite.runPart1(S, 'run_06_full_coag_turb_disagg_rate10');
DiagnosticsSuite.runPart2(S, 'run_06_full_coag_turb_disagg_rate10');

disp('RUN 06 DONE: full (coag + turb + disagg RATE)');
