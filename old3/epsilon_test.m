%% === epsilon_sweep_test.m ===
% Test the effect of turbulence (ε) on particle aggregation & export
% Runs CoagulationSimulation for multiple epsilon values and compares flux vs time
% and plots final PSD for each ε
% Author: Turjo | UGA | Nov 2025
% --------------------------------------------------------------

clear; clc; close all;

% --- ε values to test ---
eps_values = [1e-7, 1e-6, 1e-5];
colors = lines(numel(eps_values));

% --- Storage for results ---
flux_all = {};
time_all = {};
Y_final = {};
grid_all = {};

fprintf('\n=== ε sweep test (1e-7 to 1e-5) ===\n');

for i = 1:numel(eps_values)
    eps_now = eps_values(i);
    fprintf('\n--- Running case %d/%d: ε = %.1e ---\n', i, numel(eps_values), eps_now);

    % === Config ===
    cfg = SimulationConfig();
    cfg.epsilon     = eps_now;
    cfg.t_final     = 30;
    cfg.use_NPP     = true;
    cfg.NPP_profile = 'sine';
    cfg.NPP_rate    = 5e-4;

    % === Run simulation ===
    sim = CoagulationSimulation(cfg);
    out = sim.run();

    % === Compute total flux ===
    grid = DerivedGrid(cfg);
    [t, flux] = OutputGenerator.simpleFluxTimeSeries(out, grid);

    time_all{i} = t;
    flux_all{i} = flux;

    % === Save final PSD data ===
    Y_final{i} = out.concentrations(end, :)';
    grid_all{i} = grid;
end

%% === FIGURE 1: Flux vs Time (ε sweep) ===
figure('Color','w'); hold on;
for i = 1:numel(eps_values)
    plot(time_all{i}, flux_all{i}, 'LineWidth', 2, 'Color', colors(i,:));
end
xlabel('Time (days)');
ylabel('Total flux (a.u.)');
title('Effect of Turbulent Dissipation Rate (ε) on Export Flux');
legend(arrayfun(@(x) sprintf('\\epsilon = %.0e W kg^{-1}', x), eps_values, 'UniformOutput', false), ...
       'Location', 'best');
grid on;
box on;

outdir = fullfile(pwd, 'results');
if ~exist(outdir, 'dir'); mkdir(outdir); end
exportgraphics(gcf, fullfile(outdir, 'epsilon_sweep_flux.png'), 'Resolution', 300);


%% === FIGURE 2: Final PSD (size spectrum) ===
figure('Color','w'); hold on;

for i = 1:numel(eps_values)
    grid = grid_all{i};
    diam = 2 * grid.getFractalRadii() * 1e4;  % cm → µm
    loglog(diam, Y_final{i}, 'LineWidth', 2, 'Color', colors(i,:));
end

xlabel('Particle diameter (µm)');
ylabel('Concentration (a.u.)');
title('Final Particle Size Distribution vs ε');
legend(arrayfun(@(x) sprintf('\\epsilon = %.0e W kg^{-1}', x), eps_values, 'UniformOutput', false), ...
       'Location','southwest');
grid on; box on;

exportgraphics(gcf, fullfile(outdir, 'epsilon_sweep_PSD.png'), 'Resolution', 300);

fprintf('\n✅ ε-sweep complete. Figures saved to: %s\n', outdir);