% ============================================================
% run_phase1_BASELINE_coagON_disaggOFF_PPgrowth.m
% Baseline: coag=ON, disagg=OFF, PP handled by growth_mode='pp'
% ============================================================

clear; close all; clc;

root = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/26 Jan';
assert(isfolder(root), 'Root directory does not exist: %s', root);

cd(root);
addpath(genpath(root));

cfg = SimulationConfig();

% Column
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

% Bookkeeping (keep these)
cfg.state_is_biovolume = true;
cfg.export_weight      = "ones";

% PP via growth_mode='pp'
cfg.enable_pp   = false;     % keep OFF (not explicit injection)
cfg.growth_mode = 'pp';
cfg.growth      = 0.15;      % your baseline value

% Optional: make solver a bit more robust
cfg.ode_options = odeset('RelTol',1e-6,'AbsTol',1e-12,'MaxStep',0.25);

% Run diagnostics suite
[outdir, dbg] = run_Diagnostics_PHASE1(cfg);

% Summary text file
writeSummary_PHASE1(outdir, cfg, dbg);

fprintf('\nPHASE-1 done.\nFolder:\n  %s\n', outdir);

% ============================================================
function writeSummary_PHASE1(outdir, cfg, dbg)

fn  = fullfile(outdir, 'SUMMARY_PHASE1.txt');
fid = fopen(fn,'w');

fprintf(fid, 'PHASE 1 DIAGNOSTICS\n');
fprintf(fid, 'Folder: %s\n\n', outdir);

fprintf(fid, 'Config:\n');
fprintf(fid, '  use_column          = %d\n', cfg.use_column);
fprintf(fid, '  z_max               = %.2f m\n', cfg.z_max);
fprintf(fid, '  dz                  = %.2f m\n', cfg.dz);
fprintf(fid, '  t_final             = %.2f d\n', cfg.t_final);
fprintf(fid, '  enable_coag         = %d\n', cfg.enable_coag);
fprintf(fid, '  enable_disagg       = %d\n', cfg.enable_disagg);
fprintf(fid, '  state_is_biovolume  = %d\n', cfg.state_is_biovolume);
fprintf(fid, '  export_weight       = %s\n', string(cfg.export_weight));
fprintf(fid, '  enable_pp           = %d\n', cfg.enable_pp);
fprintf(fid, '  growth_mode         = %s\n', string(cfg.growth_mode));
fprintf(fid, '  growth              = %.4g\n', cfg.growth);

fprintf(fid, '\nKey checks:\n');

% Prefer the RHS-based residual if present (best)
R = [];
tag = '';
if isfield(dbg,'RES_simple_rhs') && ~isempty(dbg.RES_simple_rhs)
    R = dbg.RES_simple_rhs(:); tag = 'SIMPLE residual (rhs)';
elseif isfield(dbg,'RES_simple') && ~isempty(dbg.RES_simple)
    R = dbg.RES_simple(:); tag = 'SIMPLE residual (legacy)';
end

if isempty(R)
    fprintf(fid, '  SIMPLE residual: not found in dbg\n');
else
    ok = isfinite(R);
    if isfield(dbg,'it_valid') && ~isempty(dbg.it_valid)
        itv = dbg.it_valid(:);
        n = min(numel(R), numel(itv));
        R  = R(1:n);
        ok = ok(1:n) & (itv(1:n) > 1);
    end
    if any(ok)
        fprintf(fid, '  %s (median, max) = %.3e , %.3e\n', tag, ...
            median(abs(R(ok))), max(abs(R(ok))));
    else
        fprintf(fid, '  %s: all invalid\n', tag);
    end
end

fprintf(fid, '\nFigures:\n');
fprintf(fid, '  diagnostics_mass_balance.png\n');
fprintf(fid, '  diagnostics_export_sizefractions_500_2000.png\n');
fprintf(fid, '  diagnostics_sizespectra_depth10.png\n');
fprintf(fid, '  diagnostics_sizespectra_depth50.png\n');
fprintf(fid, '  diagnostics_sizespectra_deep.png\n');

fclose(fid);
end