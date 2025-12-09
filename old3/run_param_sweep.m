%% run_param_sweep.m  — quick sensitivity
clear; close all;

matfile = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code/epsilon_daily.mat';
outroot = '/Users/turjo/Desktop/sweeps_obs_eps';
if ~exist(outroot,'dir'), mkdir(outroot); end

eps_scales = [0.5 1 2];
fr_dims    = [2.2 2.5 2.7];
use_disagg = [false true];

[t_days, eps_series] = load_epsilon_series(matfile,0, ...
    'time_field','mtime','eps_field','eps','depth_field','z','mld_field','mld', ...
    'unit','linear','agg','mld','topN',3,'resample_daily',true,'smooth_win',3);

run_id = 0;
for s = eps_scales
  for fd = fr_dims
    for df = use_disagg
      run_id = run_id+1;
      tag   = sprintf('s%g_fd%.1f_disagg%d', s, fd, df);
      rundir = fullfile(outroot, tag); if ~exist(rundir,'dir'), mkdir(rundir); end

      cfg = SimulationConfig( ...
            'epsilon_profile','observed', ...
            'epsilon_time',t_days, ...
            'epsilon_series', s*eps_series, ...
            'epsilon_ref',1e-6, ...
            'fr_dim',fd, ...
            'disagg_use_nonlinear', df);

      sim = CoagulationSimulation(cfg);
      res = sim.run();

      % Precompute outputs for metrics/extra plots
      od = OutputGenerator.spectraAndFluxes(res.time, res.concentrations, sim.grid, sim.config);

      % ---- quick metrics for table/plots across runs ----
      metrics = table;
      metrics.run      = string(tag);
      metrics.epsScale = s;
      metrics.frDim    = fd;
      metrics.disagg   = df;
      metrics.mass0    = sum(res.concentrations(1,:));
      metrics.massEnd  = sum(res.concentrations(end,:));
      metrics.fluxMax  = max(od.total_flux);
      [~,imax]         = max(od.total_flux);
      metrics.tFluxMax = res.time(imax);
      metrics.vbar0    = od.total_flux(1) / max(od.total_mass(1), eps) / 1e6;
      metrics.vbarEnd  = od.total_flux(end)/ max(od.total_mass(end),eps)/1e6;

      T = table(res.time, od.total_flux, od.total_mass, ...
                od.total_flux ./ max(od.total_mass,eps) / 1e6, ...
                'VariableNames',{'time','total_flux','total_mass','vbar'});
      writetable(metrics, fullfile(rundir,'metrics.csv'));
      writetable(T,       fullfile(rundir,'series.csv'));

      % -------- core figures (1–4) --------
      sim.generateOutputs(true);

      % -------- bonus figures (5–6) so they always exist --------
      OutputGenerator.plotConcentrationEvolution(res.time, res.concentrations, od.total_mass); % Fig 5
      OutputGenerator.plotFluxSurface(res, sim.grid);                                          % Fig 6

      % -------- save ALL open figures (PNG, PDF, FIG) --------
      save_all_figs(rundir);

      close all
    end
  end
end

disp('✅ Sweep finished. Check folders under:'); disp(outroot);

%% ---------- helpers ----------
function save_all_figs(outdir)
% Save every open figure as PNG+PDF+.fig with stable names.
drawnow; pause(0.01);                                % ensure graphics complete
figs = sort(findall(groot,'Type','figure'), 'ascend', 'ComparisonMethod','real');

% Mapping fig number → friendly name (fallback to figN if not mapped)
labels = containers.Map( ...
    {1,                 2,                     3,                        4,                        5,                          6}, ...
    {'01_spectra_panel','02_mass_balance',     '03_coag_vs_sett_ratio',  '04_coag_sett_surface',   '05_sectional_conc_ts',     '06_flux_size_time'} );

for f = figs(:)'
    % Hide axes toolbars safely (prevents export warning)
    ax = findall(f,'Type','axes');
    for a = ax'
        if isprop(a,'Toolbar') && ~isempty(a.Toolbar)
            a.Toolbar.Visible = 'off';
        end
    end

    % Choose filename
    num = f.Number;
    if isKey(labels,num), base = labels(num); else, base = sprintf('fig%d',num); end
    pngPath = fullfile(outdir, [base '.png']);
    pdfPath = fullfile(outdir, [base '.pdf']);
    figPath = fullfile(outdir, [base '.fig']);

    % Export
    try
        exportgraphics(f, pngPath, 'Resolution', 300, 'BackgroundColor','white');
        exportgraphics(f, pdfPath, 'ContentType','vector', 'BackgroundColor','white');
        savefig(f, figPath);
    catch ME
        warning('Save failed for figure %d: %s', num, ME.message);
    end
end
end