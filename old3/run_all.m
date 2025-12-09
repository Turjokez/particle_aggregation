%% run_all.m — one-button baseline + extras
clear; close all;

matfile = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code/epsilon_daily.mat';
[t_days, eps_series] = load_epsilon_series(matfile,0, ...
    'time_field','mtime','eps_field','eps','depth_field','z','mld_field','mld', ...
    'unit','linear','agg','mld','topN',3,'resample_daily',true,'smooth_win',3);

cfg = SimulationConfig('epsilon_profile','observed','epsilon_time',t_days, ...
       'epsilon_series',eps_series,'epsilon_ref',1e-6,'disagg_use_nonlinear',true);

sim = CoagulationSimulation(cfg);
res = sim.run();
sim.generateOutputs(true);                     % Figs 1–4

od = OutputGenerator.spectraAndFluxes(res.time, res.concentrations, sim.grid, sim.config);
OutputGenerator.plotConcentrationEvolution(res.time, res.concentrations, od.total_mass);   % Fig 5
OutputGenerator.plotFluxSurface(res, sim.grid);                                            % Fig 6
% NEW (only plot if matrices exist in res)
if isfield(res,'diagnostics') && isfield(res.diagnostics,'gains') && isfield(res.diagnostics,'losses')
    od = OutputGenerator.spectraAndFluxes(res.time, res.concentrations, sim.grid, sim.config);
    OutputGenerator.plotGainsLossesSurface(res, od);
end

outdir = '/Users/turjo/Desktop/run_obs_eps';
od = OutputGenerator.spectraAndFluxes(res.time, res.concentrations, sim.grid, sim.config);
OutputGenerator.plotConcentrationEvolution(res.time, res.concentrations, od.total_mass);   % Fig 5
OutputGenerator.plotFluxSurface(res, sim.grid);                                           % Fig 6

if ~exist(outdir,'dir'), mkdir(outdir); end
arrayfun(@(k) exportgraphics(k, fullfile(outdir, sprintf('fig%d.png',k.Number)), 'Resolution',300), findall(0,'Type','figure'));
writematrix([res.time, od.total_flux, od.total_mass], fullfile(outdir,'series_baseline.csv'));
disp(['✅ Done. Outputs saved in: ' outdir]);