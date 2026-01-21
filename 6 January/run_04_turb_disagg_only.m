% run_04_coag_turb_no_disagg.m
% Baseline: coag ON + turbulence forcing ON, disagg OFF

clear; close all;

cfg = SimulationConfig();
cfg.use_column  = true;
cfg.z_max       = 65;
cfg.dz          = 5;
cfg.n_sections  = 24;

% ---- physics ----
cfg.enable_coag   = true;
cfg.enable_disagg = false;          % IMPORTANT

cfg.growth_mode = 'shift';

% ---- turbulence forcing (pulse) ----
eps0 = 1e-8;
eps1 = 1e-6;
t1   = 10;
t2   = 12;

cfg.eps_fun = @(t,z) (eps0 + (eps1-eps0) * double(t >= t1 & t <= t2));
cfg.eps_ref = eps0;

cfg.debug_rhs_units = false;

S = CoagulationSimulation(cfg);
out = S.run();

DiagnosticsSuite.runPart1(S, 'diag_04_coag_turb_no_disagg_part1');
DiagnosticsSuite.runPart2(S, 'diag_04_coag_turb_no_disagg_part2');

disp('RUN 04 DONE: coag + turb, disagg OFF');