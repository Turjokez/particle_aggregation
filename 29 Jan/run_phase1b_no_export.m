% ============================================================
% run_phase1B_no_export.m
% PHASE-1B driver: coag=ON, disagg=OFF, PP=ON, Export=OFF
% Goal: no loss term => M(t) should increase only due to PP input
% ============================================================

clear; close all; clc;

root = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/26 Jan';
assert(isfolder(root), 'Root directory does not exist: %s', root);

cd(root);
addpath(genpath(root));

% -------------------------
% CONFIG
% -------------------------
cfg = SimulationConfig();

% Column settings
cfg.use_column = true;
cfg.z_max      = 150;
cfg.dz         = 5;

% Time
cfg.t_init  = 0;
cfg.t_final = 30;
cfg.delta_t = 1;

% Physics
cfg.enable_coag   = true;
cfg.enable_disagg = false;

% --- bookkeeping ---
cfg.state_is_biovolume = true;
cfg.export_weight      = "ones";

% --- PP explicit injection ---
cfg.growth      = 0;
cfg.growth_mode = 'shift';

cfg.enable_pp   = true;
cfg.pp_rate     = 0.1;   % [cm^3/cm^2/day] because state_is_biovolume=true
cfg.pp_bin      = 1;
cfg.pp_layer    = 1;

% -------------------------
% EXPORT OFF (sinking OFF)
% -------------------------
% One of these names should match your SimulationConfig.
% Keep them all (no deletion). Ones that don't exist will be skipped.

if isprop(cfg,'enable_sinking'); cfg.enable_sinking = false; end
if isprop(cfg,'enable_settling'); cfg.enable_settling = false; end
if isprop(cfg,'enable_export'); cfg.enable_export = false; end

% Sometimes sinking speed is set as w_sink / wsink / w
if isprop(cfg,'w_sink'); cfg.w_sink = 0; end
if isprop(cfg,'wsink');  cfg.wsink  = 0; end
if isprop(cfg,'w');      cfg.w      = 0; end

% -------------------------
% RUN DIAGNOSTICS
% -------------------------
[outdir, dbg] = run_Diagnostics_PHASE1(cfg);

% Optional summary (same local-function rule as Phase-1)
writeSummary_PHASE1(outdir, cfg, dbg);

fprintf('\nPHASE-1B done.\nFolder:\n  %s\n', outdir);

% ============================================================
% Helper: write summary (keep inside same file)
% ============================================================
function writeSummary_PHASE1(outdir, cfg, dbg)
fn = fullfile(outdir, 'SUMMARY_PHASE1B_NO_EXPORT.txt');
fid = fopen(fn,'w');

fprintf(fid, 'PHASE 1B: EXPORT OFF SANITY CHECK\n');
fprintf(fid, 'Folder: %s\n\n', outdir);

fprintf(fid, 'Config:\n');
fprintf(fid, '  use_column   = %d\n', cfg.use_column);
fprintf(fid, '  z_max        = %.2f m\n', cfg.z_max);
fprintf(fid, '  dz           = %.2f m\n', cfg.dz);
fprintf(fid, '  t_final      = %.2f d\n', cfg.t_final);
fprintf(fid, '  enable_coag  = %d\n', cfg.enable_coag);
fprintf(fid, '  enable_disagg= %d\n', cfg.enable_disagg);

if isprop(cfg,'enable_pp'),     fprintf(fid, '  enable_pp    = %d\n', cfg.enable_pp); end
if isprop(cfg,'pp_rate'),       fprintf(fid, '  pp_rate      = %.4g\n', cfg.pp_rate); end

% record export-off knobs if present
if isprop(cfg,'enable_sinking'), fprintf(fid,'  enable_sinking = %d\n', cfg.enable_sinking); end
if isprop(cfg,'w_sink'),         fprintf(fid,'  w_sink         = %.4g\n', cfg.w_sink); end
if isprop(cfg,'wsink'),          fprintf(fid,'  wsink          = %.4g\n', cfg.wsink); end
if isprop(cfg,'w'),              fprintf(fid,'  w              = %.4g\n', cfg.w); end

fprintf(fid, '\nKey expectation:\n');
fprintf(fid, '  Export rate ~ 0 for all t\n');
fprintf(fid, '  M(t) should grow ~ linearly with slope ~ pp_rate\n');
fprintf(fid, '  Residual dM/dt - (PP - Export) should be near 0\n');

fclose(fid);
end