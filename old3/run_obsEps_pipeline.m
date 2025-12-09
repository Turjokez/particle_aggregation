%% run_obsEps_pipeline.m
% One-click: load ε(t) from epsilon_daily.mat, run simulation, save figs, export CSV

clear; close all; clc;

% --- paths ----------------------------------------------------------------
matfile = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code/epsilon_daily.mat';
outdir  = '/Users/turjo/Desktop/run_obs_eps';
if ~exist(outdir,'dir'), mkdir(outdir); end

% --- 1) Load ε(t) (ML-mean, daily, linear units) --------------------------
[t_days, eps_series] = load_epsilon_series(matfile, 0, ...
    'time_field','mtime', 'eps_field','eps', ...
    'depth_field','z',    'mld_field','mld', ...
    'unit','linear', 'agg','mld', 'topN',3, ...
    'resample_daily',true, 'smooth_win',3);

fprintf('ε range: %.2e to %.2e W kg^-1\n', min(eps_series), max(eps_series));

% --- 2) Configure + run ----------------------------------------------------
cfg = SimulationConfig('epsilon_profile','observed', ...
    'epsilon_time',t_days, 'epsilon_series',eps_series, ...
    'epsilon_ref',1e-6, 'disagg_use_nonlinear',true);

sim = CoagulationSimulation(cfg);
res = sim.run();
sim.generateOutputs(true);           % Figures 1–4

% --- 3) Save figures -------------------------------------------------------
print(1, fullfile(outdir,'Fig1_spectra_flux.png'),        '-dpng','-r300');
print(2, fullfile(outdir,'Fig2_total_mass_balance.png'),  '-dpng','-r300');
print(3, fullfile(outdir,'Fig3_coag_over_sett_time.png'), '-dpng','-r300');
print(4, fullfile(outdir,'Fig4_coag_over_sett_surface.png'), '-dpng','-r300');

% --- 4) Export key time series --------------------------------------------
% Grab a grid object for export calcs
if isprop(sim,'derivedGrid')
    gridobj = sim.derivedGrid;
elseif isprop(sim,'grid')
    gridobj = sim.grid;
elseif ismethod(sim,'getGrid')
    gridobj = sim.getGrid();
else
    % rebuild the same grid from the sim's config (works with your classes)
    gridobj = DerivedGrid(sim.config);
end

% Now compute spectra/fluxes for saving
od = OutputGenerator.spectraAndFluxes(res.time, res.concentrations, gridobj, sim.config);
L = res.diagnostics.total_losses;        % <- new location
% settle field can be 'sett' or 'settl' depending on builder
if     isfield(L,'sett'),  Lsett = L.sett;
elseif isfield(L,'settl'), Lsett = L.settl;
else,  error('No settling loss field found (sett/settl).');
end
ratio_ts = L.coag ./ max(Lsett, eps);    % guarded division
ratio_ts(~isfinite(ratio_ts)) = NaN;

tbl = table(res.time(:), ratio_ts(:), od.total_flux(:), od.total_mass(:), ...
    'VariableNames',{'time_d','coag_over_sett','total_flux','total_mass'});
writetable(tbl, fullfile(outdir,'time_series_diagnostics.csv'));

% --- 5) Shade storm windows on Fig 3 (optional) ---------------------------
storms = [7 10; 14 15; 18 20; 21 22];   % edit if needed
figure(3); hold on; yl = ylim;
for k = 1:size(storms,1)
    x = storms(k,:);
    patch([x(1) x(2) x(2) x(1)], [yl(1) yl(1) yl(2) yl(2)], ...
          [0.85 0.9 1.0], 'EdgeColor','none', 'FaceAlpha',0.25);
    text(mean(x), yl(2)*0.95, sprintf('Storm %d',k), ...
         'HorizontalAlignment','center','FontSize',10);
end
hold off;
print(3, fullfile(outdir,'Fig3_with_storms.png'), '-dpng','-r300');

% --- 6) Soft-cap Fig 4 spike for nicer viewing (optional) -----------------
figure(4);
h = findobj(gca,'Type','Surface');
if ~isempty(h)
    Z = get(h,'ZData'); Zc = Z; Zc(~isfinite(Zc)) = NaN;
    hi = prctile(Zc(:), 99);      % soft cap at 99th percentile
    zlim([0 hi]); caxis([0 hi]);
    print(4, fullfile(outdir,'Fig4_surface_softcapped.png'), '-dpng','-r300');
end

% --- 7) Quick ε-sensitivity runs (optional) -------------------------------
sensdir = fullfile(outdir,'sensitivity'); if ~exist(sensdir,'dir'), mkdir(sensdir); end
base_e  = cfg.epsilon_series;
runs = { 'const',  median(base_e)*ones(size(base_e)); ...
         'low',    base_e/3; ...
         'high',   base_e*3 };
for i = 1:size(runs,1)
    tag = runs{i,1}; es = runs{i,2};
    cfg2 = cfg; cfg2.epsilon_series = es;
    s2 = CoagulationSimulation(cfg2);
    r2 = s2.run(); %#ok<NASGU>  % not used here, just for figs
    s2.generateOutputs(true);
    print(2, fullfile(sensdir, sprintf('Fig2_massbalance_%s.png',tag)), '-dpng','-r300');
    print(3, fullfile(sensdir, sprintf('Fig3_coag_over_sett_%s.png',tag)), '-dpng','-r300');
end

disp('✅ Done. Figures and CSV saved in:'); disp(outdir);