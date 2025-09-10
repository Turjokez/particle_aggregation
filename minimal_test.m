% Minimal test script for OOP coagulation simulation
fprintf('Testing basic simulation creation...\n');

try
    % Create simulation with minimal parameters
    config = SimulationConfig('n_sections', 5, 't_final', 5.0);
    fprintf('✓ Configuration created\n');
    
    % Create simulation
    sim = CoagulationSimulation(config);
    fprintf('✓ Simulation created\n');
    
    % Try to run (this might take a while)
    fprintf('Running simulation...\n');
    result = sim.run();
    fprintf('✓ Simulation completed successfully\n');
    
catch ME
    fprintf('✗ Error: %s\n', ME.message);
    fprintf('Location: %s (line %d)\n', ME.stack(1).name, ME.stack(1).line);
end
