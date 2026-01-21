% run_02_coag_only.m
clear; close all;

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 65;
cfg.dz         = 5;
cfg.n_sections = 24;

cfg.enable_coag   = true;
cfg.enable_disagg = false;

% weaken coagulation
if isprop(cfg,'coag_scale');   cfg.coag_scale   = 0.2; end
if isprop(cfg,'beta_scale');   cfg.beta_scale   = 0.2; end
if isprop(cfg,'kernel_scale'); cfg.kernel_scale = 0.2; end

S = CoagulationSimulation(cfg);
out = S.run();

DiagnosticsSuite.runPart1(S, 'diag_02_coag_only_part1');
DiagnosticsSuite.runPart2(S, 'diag_02_coag_only_part2');

disp('RUN 02 DONE: coag only (weak)');
