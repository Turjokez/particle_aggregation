% run_01_control_noCoag.m
clear; close all;

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 65;
cfg.dz         = 5;
cfg.n_sections = 24;

cfg.enable_coag   = false;
cfg.enable_disagg = false;
% no turbulence provided

S = CoagulationSimulation(cfg);
out = S.run();

DiagnosticsSuite.runPart1(S, 'diag_01_control_part1');
DiagnosticsSuite.runPart2(S, 'diag_01_control_part2');

disp('RUN 01 DONE: control (no coag, no turbulence, no disagg)');
