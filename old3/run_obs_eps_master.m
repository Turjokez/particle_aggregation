%% ========================================================================
% run_obs_eps_master.m
% Master driver: baseline + full figure + (optional) sweep + clean saving
% Author: Turjo
% ========================================================================

clear; close all; clc;

%% ------------------------------ USER PATHS ------------------------------
matfile = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code/epsilon_daily.mat';
outroot = '/Users/turjo/Desktop/run_obs_eps_master';
if ~exist(outroot,'dir'), mkdir(outroot); end

%% ------------------------------ CONTROLS --------------------------------
doBaseline = true;    % run a single baseline and export figures
doSweep    = true;    % run quick sensitivity sweep and export figures

% Baseline knobs
baseline.fr_dim   = 2.5;
baseline.disagg   = true;
baseline.epsScale = 1.0;
baseline.tag      = 'baseline';

% Sweep knobs (only used if doSweep == true)
eps_scales = [0.5 1 2];
fr_dims    = [2.2 2.5 2.7];
use_disagg = [false true];

%% ----------------------- Load observed ε(t) series ----------------------
[t_days, eps_series] = load_epsilon_series( ...
    matfile, 0, ...
    'time_field','mtime', 'eps_field','eps', ...
    'depth_field','z', 'mld_field','mld', ...
    'unit','linear', 'agg','mld', 'topN',3, ...
    'resample_daily',true, 'smooth_win',3);

%% ============================== BASELINE ================================
if doBaseline
    tag    = baseline.tag;
    rundir = fullfile(outroot, tag); if ~exist(rundir,'dir'), mkdir(rundir); end

    [res, sim, od] = run_and_save( ...
        t_days, baseline.epsScale*eps_series, baseline.fr_dim, baseline.disagg, ...
        rundir, tag);

    % Leave variables in base for quick ad-hoc plotting if desired
    assignin('base','res',res);
    assignin('base','sim',sim);
    assignin('base','od', od);
end

%% =============================== SWEEP =================================
if doSweep
    for s = eps_scales
        for fd = fr_dims
            for dg = use_disagg
                tag    = sprintf('s%g_fd%.1f_disagg%d', s, fd, dg);
                rundir = fullfile(outroot, tag); if ~exist(rundir,'dir'), mkdir(rundir); end

                [res, sim, od] = run_and_save( ...
                    t_days, s*eps_series, fd, dg, rundir, tag);

                % Quick per-run metrics table
                metrics = table;
                metrics.run      = string(tag);
                metrics.epsScale = s;
                metrics.frDim    = fd;
                metrics.disagg   = dg;
                metrics.mass0    = od.total_mass(1);
                metrics.massEnd  = od.total_mass(end);
                [mx, imx]        = max(od.total_flux);
                metrics.fluxMax  = mx;
                metrics.tFluxMax = res.time(imx);
                metrics.vbar0    = od.total_flux(1)  ./ max(od.total_mass(1),eps)  / 1e6;
                metrics.vbarEnd  = od.total_flux(end) ./ max(od.total_mass(end),eps) / 1e6;

                writetable(metrics, fullfile(rundir,'metrics.csv'));
            end
        end
    end
end

disp('✅ All done. Check outputs under:'); disp(outroot);


%% =======================================================================
%                              LOCAL FUNCTIONS
% =======================================================================

function [res, sim, od] = run_and_save(t_days, eps_series, fr_dim, use_disagg, outdir, tag)
%RUN_AND_SAVE  Build config, run, compute derived outputs, and save figures/CSVs

    cfg = make_cfg(t_days, eps_series, fr_dim, use_disagg);
    sim = CoagulationSimulation(cfg);
    res = sim.run();
    od  = OutputGenerator.spectraAndFluxes(res.time, res.concentrations, sim.grid, sim.config);

    % 1) Generate the model's standard panels and save all open figures
    sim.generateOutputs(true);
    save_all_figs(outdir, 'panel');

    % 2) Make the compact 2×3 "full figure" and save
    plot_full_figure(res, od, sim.grid, sim.config, outdir, tag);

    % 3) Time-series CSV (helpful for overlays later)
    T = table(res.time(:), od.total_flux(:), od.total_mass(:), ...
              od.total_flux(:)./max(od.total_mass(:),eps)/1e6, ...
              'VariableNames',{'time','total_flux','total_mass','vbar_m_per_d'});
    writetable(T, fullfile(outdir,'series.csv'));

    % 4) Optional: try 3-D diagnostics if available (safe wrapper)
    try
        % If your codebase has this method, it will produce a 3-D surface.
        OutputGenerator.plotGainsLossesSurface(res);
        save_all_figs(outdir, 'diag');
    catch
        % Quietly skip when diagnostics aren't present
    end

    close all;
end


function cfg = make_cfg(t_days, eps_series, fr_dim, use_disagg)
%MAKE_CFG  Small convenience builder for SimulationConfig

    cfg = SimulationConfig( ...
        'epsilon_profile',       'observed', ...
        'epsilon_time',          t_days, ...
        'epsilon_series',        eps_series, ...
        'epsilon_ref',           1e-6, ...
        'fr_dim',                fr_dim, ...
        'disagg_use_nonlinear',  logical(use_disagg) ...
    );
end


function plot_full_figure(res, od, gridObj, cfg, outdir, tag)
%PLOT_FULL_FIGURE  2×3 dashboard: spectra, conc(t), flux spectra, flux(t), <v>(t)
% Robust to different field names in 'od' and falls back to grid if needed.

    if ~exist(outdir,'dir'), mkdir(outdir); end

    % ---- tolerant getters ------------------------------------------------
    t            = res.time(:);
    Dcm          = get_diam_cm(od, gridObj);                             % <-- fix
    n_init       = get_field(od, {'n_init','n0','num_init','number_init','n_initial'});
    conc_t       = get_field(od, {'sectional_conc','conc_sectional','sectionalConcentration','sect_conc'});
    flux_spec    = get_field(od, {'flux_spectra','flux_by_size','flux_spectrum','fluxSpec'});
    flux_sect_t  = get_field(od, {'sectional_flux','flux_sectional','sect_flux'});
    total_flux   = get_field(od, {'total_flux','flux_total','FluxTotal'});
    total_mass   = get_field(od, {'total_mass','mass_total','MassTotal'});

    % safe average velocity (m d^-1)
    vbar = total_flux ./ max(total_mass, eps) / 1e6;

    % ---- figure layout ---------------------------------------------------
    fh = figure('Color','w','Position',[100 80 1300 720]);
    tl = tiledlayout(fh,2,3,'TileSpacing','compact','Padding','compact');

    % (1) Number spectrum (initial vs D^-4)
    nexttile(tl,1);
    if ~isempty(Dcm) && ~isempty(n_init)
        loglog(Dcm, n_init, 'r','LineWidth',1.5); hold on;
        ref = (Dcm.^-4); ref = ref * (n_init(1)/max(ref));
        loglog(Dcm, ref, 'b');
        xlabel('Particle diameter [cm]');
        ylabel('Number spectrum [# cm^{-4}]');
        grid on;
    else
        text(0.5,0.5,'(missing n\_init or D)','HorizontalAlignment','center');
        axis off;
    end

    % (2) Sectional concentration vs time (+ total mass)
    nexttile(tl,2);
    if ~isempty(conc_t)
        plot(t, conc_t,'LineWidth',1); set(gca,'YScale','log'); hold on;
        if ~isempty(total_mass), plot(t, total_mass,'m:*','LineWidth',1.4); end
        xlabel('Time [d]'); ylabel('Sectional concentration [vol/vol/sect]');
        grid on;
    else
        text(0.5,0.5,'(missing sectional conc)','HorizontalAlignment','center'); axis off;
    end

    % (3) Volume flux spectra (by size)
    nexttile(tl,4);
    if ~isempty(Dcm) && ~isempty(flux_spec)
        plot(Dcm, flux_spec,'LineWidth',1);
        xlabel('Particle image diameter [cm]');
        ylabel('Volume flux spectra [cm^2 m^{-2} d^{-1}]');
        grid on;
    else
        text(0.5,0.5,'(missing flux spectra or D)','HorizontalAlignment','center'); axis off;
    end

    % (4) Sectional flux time series (+ total)
    nexttile(tl,5);
    if ~isempty(flux_sect_t)
        plot(t, flux_sect_t,'LineWidth',1); hold on;
        if ~isempty(total_flux), plot(t, total_flux,'m:*','LineWidth',1.6); end
        xlabel('Time [d]'); ylabel('Sectional Flux [cm^3 m^{-2} d^{-1} sect^{-1}]');
        grid on;
    else
        text(0.5,0.5,'(missing sectional flux)','HorizontalAlignment','center'); axis off;
    end

    % (5) Average settling velocity
    nexttile(tl,6);
    if ~isempty(total_flux) && ~isempty(total_mass)
        plot(t, vbar, 'LineWidth',1.6);
        xlabel('Time [d]'); ylabel('Average v [m d^{-1}]'); grid on;
    else
        text(0.5,0.5,'(missing total flux/mass)','HorizontalAlignment','center'); axis off;
    end

    title(tl, sprintf('Full figure — %s', tag), 'FontWeight','bold');

    hide_axes_toolbars(fh);
    exportgraphics(fh, fullfile(outdir, sprintf('full_%s.png', tag)), 'Resolution', 300);
end

% ---------------------------- helpers ------------------------------------

function val = get_field(S, candidates)
% Return the first present field in struct/object S from 'candidates'.
    val = [];
    for k = 1:numel(candidates)
        name = candidates{k};
        try
            if isstruct(S) && isfield(S,name)
                val = S.(name);  return;
            elseif isobject(S) && isprop(S,name)
                val = S.(name);  return;
            end
        catch, end
    end
end

function Dcm = get_diam_cm(od, gridObj)
% Try several names in 'od', then in 'gridObj'. As last resort, infer length.
    Dcm = get_field(od, {'diam_cm','Dcm','diameter_cm','diam','D'});
    if ~isempty(Dcm), Dcm = Dcm(:); return; end

    % look on grid object
    Dcm = get_field(gridObj, {'diam_cm','Dcm','diameter_cm','Dcenters_cm','Dcenters'});
    if ~isempty(Dcm), Dcm = Dcm(:); return; end

    % last resort: infer from flux spectra width if available
    Fspec = get_field(od, {'flux_spectra','flux_by_size','flux_spectrum','fluxSpec'});
    if ~isempty(Fspec)
        n = size(Fspec,2);
        Dcm = logspace(log10(1e-3), log10(1), n).';   % reasonable placeholder
    else
        Dcm = [];
    end
end

function save_all_figs(outdir, prefix)
%SAVE_ALL_FIGS  Save every open figure as PNG, ordered by figure number.

    if nargin < 2, prefix = 'fig'; end
    if ~exist(outdir,'dir'), mkdir(outdir); end

    figs = findall(groot,'Type','figure');
    if isempty(figs), return; end

    nums = arrayfun(@(f) f.Number, figs);
    [~, ord] = sort(nums);      % robust sort for numeric figure handles
    figs = figs(ord);

    for f = figs'
        hide_axes_toolbars(f);
        fn = sprintf('%s_%02d.png', prefix, f.Number);
        exportgraphics(f, fullfile(outdir, fn), 'Resolution', 300);
    end
end


function hide_axes_toolbars(fh)
%HIDE_AXES_TOOLBARS  Avoid export warning by hiding per-axes toolbars.

    ax = findall(fh,'Type','axes');
    for a = ax'
        if isprop(a,'Toolbar') && ~isempty(a.Toolbar)
            a.Toolbar.Visible = 'off';
        end
    end
end