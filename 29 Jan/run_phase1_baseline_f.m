% ============================================================
% run_phase1_baseline.m
% cfg.enable_coag = false;
% PHASE-1 ONLY driver: baseline run -> 3 diagnostics deliverables
%
% Creates a folder results_phase1_YYYYMMDD_HHMMSS with:
%   - cfg_snapshot.mat
%   - out_baseline.mat
%   - dbg.mat
%   - diagnostics_mass_balance.png
%   - diagnostics_export_sizefractions_500_2000.png
%   - diagnostics_sizespectra_depth10.png
%   - diagnostics_sizespectra_depth50.png
%   - diagnostics_sizespectra_deep.png
% ============================================================

clear; close all; clc;

root = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/26 Jan';

assert(isfolder(root), 'Root directory does not exist: %s', root);

cd(root);
addpath(genpath(root));

% -------------------------
% BASELINE CONFIG (EDIT ONLY THIS BLOCK)
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

% Physics (set this to your trusted baseline)
cfg.enable_coag   = false;
cfg.enable_disagg = false;

% Choose ONE PP approach (keep whatever matched your "good" run)
% Option A: PP as growth (recommended if that is your baseline)
cfg.enable_pp   = false;
cfg.growth_mode = 'pp';    % <-- if your baseline used this
% cfg.growth    = 0.1;     % only if your model needs it

% Option B: Explicit PP injection (only if baseline used enable_pp=true)
% cfg.growth      = 0;
% cfg.growth_mode = 'shift';
% cfg.enable_pp   = true;
% cfg.pp_rate     = 0.1;
% cfg.pp_bin      = 1;
% cfg.pp_layer    = 1;

% -------------------------
% RUN PHASE-1 DIAGNOSTICS
% -------------------------
[outdir, dbg] = run_Diagnostics_PHASE1(cfg);

% Write short summary
writeSummary_PHASE1(outdir, cfg, dbg);

fprintf('\nPHASE-1 done.\nFolder:\n  %s\n', outdir);

% ============================================================
% Helper: write summary
% ============================================================
function writeSummary_PHASE1(outdir, cfg, dbg)
fn = fullfile(outdir, 'SUMMARY_PHASE1.txt');
fid = fopen(fn,'w');

fprintf(fid, 'PHASE 1 BASELINE DIAGNOSTICS\n');
fprintf(fid, 'Folder: %s\n\n', outdir);

fprintf(fid, 'Config:\n');
fprintf(fid, '  use_column   = %d\n', cfg.use_column);
fprintf(fid, '  z_max        = %.2f m\n', cfg.z_max);
fprintf(fid, '  dz           = %.2f m\n', cfg.dz);
fprintf(fid, '  t_final      = %.2f d\n', cfg.t_final);
fprintf(fid, '  enable_coag  = %d\n', cfg.enable_coag);
fprintf(fid, '  enable_disagg= %d\n', cfg.enable_disagg);

if isprop(cfg,'enable_pp')
    fprintf(fid, '  enable_pp    = %d\n', cfg.enable_pp);
end
if isprop(cfg,'growth_mode')
    fprintf(fid, '  growth_mode  = %s\n', string(cfg.growth_mode));
end
if isprop(cfg,'growth')
    fprintf(fid, '  growth       = %.4g\n', cfg.growth);
end
if isprop(cfg,'pp_rate')
    fprintf(fid, '  pp_rate      = %.4g\n', cfg.pp_rate);
    fprintf(fid, '  pp_bin       = %d\n', cfg.pp_bin);
    fprintf(fid, '  pp_layer     = %d\n', cfg.pp_layer);
end

fprintf(fid, '\nKey checks:\n');
ok = isfinite(dbg.RES_simple) & (dbg.it_valid > 1);
fprintf(fid, '  SIMPLE residual (median, max) = %.3e , %.3e\n', ...
    median(abs(dbg.RES_simple(ok))), max(abs(dbg.RES_simple(ok))));

fprintf(fid, '\nFigures:\n');
fprintf(fid, '  diagnostics_mass_balance.png\n');
fprintf(fid, '  diagnostics_export_sizefractions_500_2000.png\n');
fprintf(fid, '  diagnostics_sizespectra_depth10.png\n');
fprintf(fid, '  diagnostics_sizespectra_depth50.png\n');
fprintf(fid, '  diagnostics_sizespectra_deep.png\n');

fclose(fid);
end