% run_coag_tests.m  (RUN 1 ONLY)

cfg = SimulationConfig();
cfg.use_column  = true;
cfg.z_max       = 65;
cfg.dz          = 5;
cfg.n_sections  = 24;

% keep disagg OFF
if isprop(cfg,'enable_disagg'); cfg.enable_disagg = false; end

% RUN 1: Coag OFF
cfg1 = cfg;
if isprop(cfg1,'enable_coag'); cfg1.enable_coag = false; end

S1   = CoagulationSimulation(cfg1);
out1 = S1.run(); %#ok<NASGU>

DiagnosticsSuite.runPart1(S1, fullfile(pwd,'diagnostics_RUN1_coagOFF_part1'));
DiagnosticsSuite.runPart2(S1, fullfile(pwd,'diagnostics_RUN1_coagOFF_part2'));

disp('Done: check Part2 coag remainder panels (should be ~0).');
