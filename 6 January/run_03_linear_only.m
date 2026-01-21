% run_03_linear_only.m
% Goal: Linear physics only (growth/sinking/transport), no coag, no disagg.

clear; close all;

cfg = SimulationConfig();
cfg.use_column  = true;
cfg.z_max       = 65;
cfg.dz          = 5;
cfg.n_sections  = 24;

% Physics switches
if isprop(cfg,'enable_coag');   cfg.enable_coag   = false; end
if isprop(cfg,'enable_disagg'); cfg.enable_disagg = false; end

% Keep your baseline choices
cfg.growth_mode = 'shift';
% cfg.growth    = 0;   % keep whatever you want; 0 means no growth injection

% Turbulence OFF (make sure nothing is feeding eps)
if isprop(cfg,'eps_fun');        cfg.eps_fun = []; end
if isprop(cfg,'epsilon_time');   cfg.epsilon_time = []; end
if isprop(cfg,'epsilon_series'); cfg.epsilon_series = []; end
if isprop(cfg,'epsilon_const');  cfg.epsilon_const = []; end

% Optional (turn OFF if you don't want huge log PDFs)
if isprop(cfg,'debug_rhs_units'); cfg.debug_rhs_units = true; end

S = CoagulationSimulation(cfg);
out = S.run();

DiagnosticsSuite.runPart1(S, 'diag_03_linear_only_part1');
DiagnosticsSuite.runPart2(S, 'diag_03_linear_only_part2');

disp('RUN 03 DONE: linear only');
