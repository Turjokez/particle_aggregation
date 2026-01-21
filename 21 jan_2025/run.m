% run.m
% One script to generate the key figures + a short summary

clear; close all; clc;

root = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/old backup14 jan';
cd(root);
addpath(genpath(root));

% RUN A: "Debug run" (PP OFF) -> show FULL mass closure + detective plots
cfg = SimulationConfig();

% Column settings 
cfg.use_column = true;
cfg.z_max      = 150;
cfg.dz         = 5;

% Time
cfg.t_init  = 0;
cfg.t_final = 30;
cfg.delta_t = 1;

% Physics toggles (keep simple)
cfg.enable_coag   = true;
cfg.enable_disagg = false;
cfg.enable_pp     = false;     % IMPORTANT: PP OFF for this run
cfg.growth_mode   = 'shift';   % ok (but not a true PP source)
% cfg.growth      = 0.15;      % keep as your default, or set to 0 if you want no growth

[outdirA, dbgA] = run_Diagnostics_FINAL(cfg);

% Save a short text summary for 
writeSummary(outdirA, cfg, dbgA, "RUN A (PP OFF)");

fprintf('\nRUN A done: %s\n', outdirA);

% RUN B: "PP test" (PP ON) -> show the PP scaling / simple-budget problem
cfg2 = SimulationConfig();

cfg2.use_column = true;
cfg2.z_max      = 150;
cfg2.dz         = 5;

cfg2.t_init  = 0;
cfg2.t_final = 30;
cfg2.delta_t = 1;

cfg2.enable_coag   = true;
cfg2.enable_disagg = false;

% Turn OFF legacy growth shift so we isolate the explicit PP source
cfg2.growth      = 0;
cfg2.growth_mode = 'shift';

% Turn ON explicit PP source
cfg2.enable_pp = true;
cfg2.pp_rate   = 0.1;     % your current value
cfg2.pp_bin    = 1;       % add in smallest bin (or choose a bin you want)
cfg2.pp_layer  = 1;       % surface layer only

[outdirB, dbgB] = run_Diagnostics_FINAL(cfg2);

writeSummary(outdirB, cfg2, dbgB, "RUN B (PP ON, growth OFF)");

fprintf('\nRUN B done: %s\n', outdirB);

% ---------------------------------------
% Quick printout in MATLAB command window
% ---------------------------------------
disp("==== QUICK CHECKS ====");
disp("RUN A: FULL residual (median, max)"); 
disp([median(abs(dbgA.RES_full(isfinite(dbgA.RES_full)))) , max(abs(dbgA.RES_full(isfinite(dbgA.RES_full))))]);

disp("RUN B: SIMPLE residual (median, max)");
disp([median(abs(dbgB.RES_simple(isfinite(dbgB.RES_simple)))) , max(abs(dbgB.RES_simple(isfinite(dbgB.RES_simple))))]);

disp("RUN B: dMdt_pp range (min, max)");
disp([min(dbgB.dMdt_pp) max(dbgB.dMdt_pp)]);



% Helper
function writeSummary(outdir, cfg, dbg, label)
    fn = fullfile(outdir, 'SUMMARY.txt');
    fid = fopen(fn, 'w');

    fprintf(fid, '%s\n', label);
    fprintf(fid, 'Folder: %s\n\n', outdir);

    fprintf(fid, 'Config:\n');
    fprintf(fid, '  use_column = %d\n', cfg.use_column);
    fprintf(fid, '  z_max      = %.2f m\n', cfg.z_max);
    fprintf(fid, '  dz         = %.2f m\n', cfg.dz);
    fprintf(fid, '  t_final    = %.2f d\n', cfg.t_final);
    fprintf(fid, '  enable_coag   = %d\n', cfg.enable_coag);
    fprintf(fid, '  enable_disagg = %d\n', cfg.enable_disagg);
    fprintf(fid, '  enable_pp     = %d\n', cfg.enable_pp);
    fprintf(fid, '  growth_mode   = %s\n', string(cfg.growth_mode));
    fprintf(fid, '  growth        = %.4g\n', cfg.growth);

    if isprop(cfg,'pp_rate')
        fprintf(fid, '  pp_rate       = %.4g\n', cfg.pp_rate);
        fprintf(fid, '  pp_bin        = %d\n', cfg.pp_bin);
        fprintf(fid, '  pp_layer      = %d\n', cfg.pp_layer);
    end

    fprintf(fid, '\nKey diagnostics:\n');
    okS = isfinite(dbg.RES_simple);
    okF = isfinite(dbg.RES_full);

    fprintf(fid, '  SIMPLE residual: median=%.3e  max=%.3e\n', ...
        median(abs(dbg.RES_simple(okS))), max(abs(dbg.RES_simple(okS))));
    fprintf(fid, '  FULL residual:   median=%.3e  max=%.3e\n', ...
        median(abs(dbg.RES_full(okF))), max(abs(dbg.RES_full(okF))));

    if isfield(dbg,'closure_max_relinf')
        fprintf(fid, '  Closure check: max rel inf=%.3e  max abs inf=%.3e\n', ...
            dbg.closure_max_relinf, dbg.closure_max_absinf);
    end

    if isfield(dbg,'EX')
        fprintf(fid, '  EX range: min=%.3e  max=%.3e\n', min(dbg.EX), max(dbg.EX));
    end
    if isfield(dbg,'dMdt_fd')
        v = dbg.dMdt_fd(isfinite(dbg.dMdt_fd));
        fprintf(fid, '  dMdt_fd range: min=%.3e  max=%.3e\n', min(v), max(v));
    end
    if isfield(dbg,'dMdt_pp')
        fprintf(fid, '  dMdt_pp range: min=%.3e  max=%.3e\n', min(dbg.dMdt_pp), max(dbg.dMdt_pp));
    end

    fprintf(fid, '\nFigures to show:\n');
    fprintf(fid, '  mass_balance_residual_full.png\n');
    fprintf(fid, '  mass_balance_residual_simple.png\n');
    fprintf(fid, '  export_sizefractions_500_1000.png\n');
    fprintf(fid, '  export_sizefractions_500_2000.png\n');
    fprintf(fid, '  fig_spectra/sizespectra_all_depths.png\n');
    fprintf(fid, '  fig_tendencies/dQdt_total_all_depths.png\n');
    fprintf(fid, '  fig_tendencies/dQdt_coag_all_depths.png\n');
    fprintf(fid, '  fig_tendencies/dQdt_linear_transport_sink_all_depths.png\n');
    fprintf(fid, '  fig_tendencies/dQdt_PP_all_depths.png (only meaningful when PP ON)\n');

    fclose(fid);
end