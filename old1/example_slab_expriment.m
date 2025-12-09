%% === example_slab_experiment.m ===
% Slab coagulation–settling experiment following Adrian's outline

% -------------------------------------------------------------
clear; clc; close all;
fprintf('\n=== Example Slab Experiment (Step 6–7, full set) ===\n');

%% === 1) CONFIGURATION ===
cfg = SimulationConfig();

cfg.t_final     = 30;       % [days]
cfg.use_NPP     = true;     
cfg.NPP_rate    = 5e-4;     
cfg.NPP_profile = 'sine';   % or 'constant'
cfg.epsilon     = 1e-6;     % [W kg^-1]
cfg.temp        = 20;       
cfg.salt        = 35;       

fprintf('Configuration ready → running simulation...\n');
sim  = CoagulationSimulation(cfg);
out  = sim.run();
grid = DerivedGrid(cfg);

fprintf('Simulation finished: %d time steps × %d sections\n',...
    numel(out.time), size(out.concentrations,2));

%% === 2) FIGURE 01 — NPP vs Time ===
figure('Color','w');
tt = linspace(0, cfg.t_final, 100);
if strcmpi(cfg.NPP_profile, 'sine')
    NPP_t = cfg.NPP_rate * (1 + 0.5*sin(2*pi*tt/cfg.t_final));
else
    NPP_t = cfg.NPP_rate * ones(size(tt));
end
plot(tt, NPP_t, 'b', 'LineWidth', 2);
xlabel('Time (days)'); ylabel('NPP (model units d^{-1})');
title('NPP vs Time'); grid on;
exportgraphics(gcf,'results/Fig_01_NPP_vs_Time.png','Resolution',300);

%% === Figure 02 — ε vs Time (log) ===
figure('Color','w');
semilogy(tt, cfg.epsilon*ones(size(tt)), 'r', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('\epsilon (W kg^{-1})');
title('Turbulent Energy Dissipation vs Time');
grid on;
ylim([1e-8 1e-5]);              % ensures axis covers full log range
yticks([1e-8 1e-7 1e-6 1e-5]);  % clean tick marks
exportgraphics(gcf,'results/Fig_02_Epsilon_vs_Time.png','Resolution',300);

%% === 4) FIGURE 03 — Flux vs Time ===
figure('Color','w');
[t, flux] = OutputGenerator.simpleFluxTimeSeries(out, grid);
title('Flux vs Time');
exportgraphics(gcf,'results/Fig_03_Flux_vs_Time.png','Resolution',300);

%% === 5) FIGURE 04 — Flux vs Size vs Time (3D surface) ===
figure('Color','w');
OutputGenerator.plotFluxSurface(out, grid);
exportgraphics(gcf,'results/Fig_04_Flux_vs_Size_vs_Time.png','Resolution',300);

%% === 6) FIGURE 05 — Multi-panel Spectra & Flux Diagnostics ===
figure('Color','w');
output_data = OutputGenerator.spectraAndFluxes(out.time, out.concentrations, grid, cfg);
OutputGenerator.plotSpectraAndFluxes(out.time, out.concentrations, output_data);
exportgraphics(gcf,'results/Fig_05_Multipanel_Spectra_and_Flux.png','Resolution',300);

