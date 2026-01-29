% ============================================================
% run_phase1_coagON_disaggOFF_PPON.m
% PHASE-1 driver: coag=ON, disagg=OFF, PP=ON (via growth_mode='pp')
% ============================================================

clear; close all; clc;

root = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/26 Jan';
assert(isfolder(root), 'Root directory does not exist: %s', root);

cd(root);
addpath(genpath(root));

% -------------------------
% CONFIG (EDIT ONLY THIS BLOCK)
% -------------------------
cfg = SimulationConfig();

% Column
cfg.use_column = true;
cfg.z_max      = 150;
cfg.dz         = 5;

% Time
cfg.t_init  = 0;
cfg.t_final = 30;
cfg.delta_t = 1;

% Physics switches
cfg.enable_coag   = true;
cfg.enable_disagg = false;

% --- critical bookkeeping (keep these) ---
cfg.state_is_biovolume = true;
cfg.export_weight      = "ones";

% -------------------------
% PP = ON (recommended way)
% -------------------------
cfg.enable_pp   = false;   % IMPORTANT: keep false for growth_mode='pp'
cfg.growth_mode = 'pp';

% If your model needs an explicit growth value, set it here:
% cfg.growth = 0.1;   % <-- only if required in your code

% -------------------------
% RUN PHASE-1 DIAGNOSTICS
% -------------------------
[outdir, dbg] = run_Diagnostics_PHASE1(cfg);

% Save short summary
writeSummary_PHASE1(outdir, cfg, dbg);

fprintf('\nPHASE-1 done.\nFolder:\n  %s\n', outdir);

% ============================================================
% Helper: write summary
% ============================================================
function writeSummary_PHASE1(outdir, cfg, dbg)
fn = fullfile(outdir, 'SUMMARY_PHASE1.txt');
fid = fopen(fn,'w');

fprintf(fid, 'PHASE 1 DIAGNOSTICS\n');
fprintf(fid, 'Folder: %s\n\n', outdir);

fprintf(fid, 'Config:\n');
fprintf(fid, '  use_column    = %d\n', cfg.use_column);
fprintf(fid, '  z_max         = %.2f m\n', cfg.z_max);
fprintf(fid, '  dz            = %.2f m\n', cfg.dz);
fprintf(fid, '  t_final       = %.2f d\n', cfg.t_final);
fprintf(fid, '  enable_coag   = %d\n', cfg.enable_coag);
fprintf(fid, '  enable_disagg = %d\n', cfg.enable_disagg);

if isprop(cfg,'enable_pp'),   fprintf(fid, '  enable_pp     = %d\n', cfg.enable_pp); end
if isprop(cfg,'growth_mode'), fprintf(fid, '  growth_mode   = %s\n', string(cfg.growth_mode)); end
if isprop(cfg,'growth'),      fprintf(fid, '  growth        = %.4g\n', cfg.growth); end

if isprop(cfg,'state_is_biovolume')
    fprintf(fid, '  state_is_biovolume = %d\n', cfg.state_is_biovolume);
end
if isprop(cfg,'export_weight')
    fprintf(fid, '  export_weight      = %s\n', string(cfg.export_weight));
end

fprintf(fid, '\nKey checks:\n');

% Pick residual that exists
if isfield(dbg,'RES_simple_rhs')
    R = dbg.RES_simple_rhs(:);
elseif isfield(dbg,'RES_simple')
    R = dbg.RES_simple(:);
else
    R = NaN;
end

ok = isfinite(R);

if isfield(dbg,'it_valid') && ~isempty(dbg.it_valid)
    itv = dbg.it_valid(:);
    n = min(numel(R), numel(itv));
    ok = ok(1:n) & (itv(1:n) > 1);
    R = R(1:n);
end

fprintf(fid, '  SIMPLE residual (median, max) = %.3e , %.3e\n', ...
    median(abs(R(ok))), max(abs(R(ok))));

fprintf(fid, '\nFigures:\n');
fprintf(fid, '  diagnostics_mass_balance.png\n');
fprintf(fid, '  diagnostics_export_sizefractions_500_2000.png\n');
fprintf(fid, '  diagnostics_sizespectra_depth10.png\n');
fprintf(fid, '  diagnostics_sizespectra_depth50.png\n');
fprintf(fid, '  diagnostics_sizespectra_deep.png\n');

fclose(fid);
end