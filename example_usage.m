%% Example Usage of OOP Coagulation Simulation
% This script demonstrates how to use the new object-oriented classes
% to run a coagulation simulation

clear all; close all;

%% Create simulation with default parameters
fprintf('Creating simulation with default parameters...\n');
sim = CoagulationSimulation();

%% UPDATED for Δx = 2 m - display vertical grid information
fprintf('Vertical grid mode: %s\n', sim.verticalMode);
if isfield(sim.verticalGrid, 'delta_z')
    fprintf(' - Surface spacing Δz = %.1f m\n', sim.verticalGrid.delta_z(1));
    fprintf(' - Number of layers: %d (max depth = %.0f m)\n', ...
        numel(sim.verticalGrid.delta_z), sim.verticalGrid.edges(end));
end

%% UPDATED for Δx = 2 m - optional stretched grid example (commented)
% sim.setVerticalMode('stretched', ...
%     'StretchedBreaks', [0 250 600 1000], ...
%     'StretchedSpacing', [25 50 100]);
% fprintf('Switched to stretched grid (%d layers).\n', numel(sim.verticalGrid.delta_z));

%% Run simulation
fprintf('Running simulation...\n');
result = sim.run();

%% Generate outputs and plots
fprintf('Generating outputs...\n');
sim.generateOutputs(true);

%% UPDATED for Δx = 2 m - report animation output
if isfield(sim.result, 'vertical_profiles') && ...
        isfield(sim.result.vertical_profiles, 'animation_file') && ...
        ~isempty(sim.result.vertical_profiles.animation_file)
    fprintf('Saved F(z) animation to: %s\n', sim.result.vertical_profiles.animation_file);
end



% UNCOMMENT EVERYTHING BELOW TO GENERATE CUSTOM CONFIG SIMULATIONS
% %% Example with custom parameters
% fprintf('\n=== Example with custom parameters ===\n');

% % Create custom configuration
% custom_config = SimulationConfig(...
%     'n_sections', 15, ...           % Fewer sections for faster computation
%     't_final', 10.0, ...            % Shorter simulation time
%     'growth', 0.1, ...              % Lower growth rate
%     'gamma', 0.05, ...              % Lower shear rate
%     'alpha', 0.8 ...                % Lower stickiness
%     );

% % Create simulation with custom config
% sim_custom = CoagulationSimulation(custom_config);

% % Run with custom time span
% tspan_custom = 0:0.5:8;  % More frequent output
% result_custom = sim_custom.run('tspan', tspan_custom);

% % Generate outputs without plots
% sim_custom.generateOutputs(false);

% %% Example of accessing individual components
% fprintf('\n=== Accessing individual components ===\n');

% % Access configuration
% fprintf('Number of sections: %d\n', sim.config.n_sections);
% fprintf('Simulation time: %.1f to %.1f days\n', sim.config.t_init, sim.config.t_final);

% % Access grid properties
% fprintf('Volume range: [%.2e, %.2e] cm³\n', ...
%     min(sim.grid.v_lower), max(sim.grid.v_upper));

% % Access results
% fprintf('Final total mass: %.2e\n', result.output_data.total_mass(end));
% fprintf('Number of time points: %d\n', length(result.time));

% %% Example of custom initial conditions
% fprintf('\n=== Custom initial conditions ===\n');

% % Create exponential initial spectrum
% v0_exp = InitialSpectrumBuilder.exponentialSpectrum(sim.config, sim.grid, 1000, 0.2);

% % Run simulation with custom initial conditions
% result_exp = sim.run('v0', v0_exp);

% fprintf('Custom initial conditions simulation completed.\n');

% %% Example of solver options
% fprintf('\n=== Custom solver options ===\n');

% % Create solver with custom options
% custom_solver = ODESolver('ode23s');
% custom_solver.setOptions('RelTol', 1e-10, 'AbsTol', 1e-20);

% % Run with custom solver
% result_custom_solver = sim.run('solver_options', custom_solver.options);

% fprintf('Custom solver simulation completed.\n');

% %% Example of analyzing specific components
% fprintf('\n=== Component analysis ===\n');

% % Analyze beta matrices
% if isfield(result, 'betas')
%     result.betas.displaySummary();
% end

% % Analyze mass balance
% if isfield(result.diagnostics, 'sectional_gains')
%     MassBalanceAnalyzer.displayBalanceSummary(...
%         result.diagnostics.sectional_gains, ...
%         result.diagnostics.sectional_losses, ...
%         result.time);
% end

% %% Example of exporting specific data
% fprintf('\n=== Data export ===\n');

% % Export just the output data
% OutputGenerator.exportData(result.output_data, 'my_output_data.mat');

% % Export full simulation object
% save('my_full_simulation.mat', 'sim');

% fprintf('Data export completed.\n');

% %% Summary
% fprintf('\n=== Summary ===\n');
% fprintf('Successfully demonstrated OOP coagulation simulation classes:\n');
% fprintf('- SimulationConfig: Parameter management\n');
% fprintf('- DerivedGrid: Precomputed constants and grids\n');
% fprintf('- KernelLibrary: Coagulation kernel implementations\n');
% fprintf('- LinearProcessBuilder: Growth, sinking, disaggregation operators\n');
% fprintf('- InitialSpectrumBuilder: Initial condition generation\n');
% fprintf('- BetaAssembler: Coagulation kernel matrix computation\n');
% fprintf('- CoagulationRHS: ODE right-hand side evaluation\n');
% fprintf('- ODESolver: Time integration\n');
% fprintf('- MassBalanceAnalyzer: Diagnostics and analysis\n');
% fprintf('- OutputGenerator: Visualization and data export\n');
% fprintf('- CoagulationSimulation: Main simulation controller\n');

% fprintf('\nThe OOP structure provides:\n');
% fprintf('- Better organization and maintainability\n');
% fprintf('- Reusable components\n');
% fprintf('- Clear separation of concerns\n');
% fprintf('- Easier testing and debugging\n');
% fprintf('- Better code documentation\n');
