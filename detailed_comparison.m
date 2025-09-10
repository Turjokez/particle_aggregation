%% Detailed Comparison of Legacy vs OOP with Plot Generation
% This script runs both versions and generates plots for visual comparison

clear all; close all;

fprintf('=== Detailed Comparison with Plot Generation ===\n\n');

%% Run Legacy Version
fprintf('Running Legacy Version...\n');
try
    % Run legacy version
    [p, opt] = SetUpCoag;

    % Calculate kernels
    b_brown = CalcBetas(p);
    b_brown.b1 = b_brown.b1*p.conBr*p.day_to_sec;
    b_brown.b2 = b_brown.b2*p.conBr*p.day_to_sec;
    b_brown.b3 = b_brown.b3*p.conBr*p.day_to_sec;
    b_brown.b4 = b_brown.b4*p.conBr*p.day_to_sec;
    b_brown.b5 = b_brown.b5*p.conBr*p.day_to_sec;

    p.kernel='KernelCurSh';
    b_shear = CalcBetas(p);
    b_shear.b1 = b_shear.b1*p.gamma*p.day_to_sec;
    b_shear.b2 = b_shear.b2*p.gamma*p.day_to_sec;
    b_shear.b3 = b_shear.b3*p.gamma*p.day_to_sec;
    b_shear.b4 = b_shear.b4*p.gamma*p.day_to_sec;
    b_shear.b5 = b_shear.b5*p.gamma*p.day_to_sec;
    b_shear.b25 = b_shear.b25*p.gamma*p.day_to_sec;

    p.kernel='KernelCurDS';
    b_ds = CalcBetas(p);
    b_ds.b1 = b_ds.b1*p.setcon*p.day_to_sec;
    b_ds.b2 = b_ds.b2*p.setcon*p.day_to_sec;
    b_ds.b3 = b_ds.b3*p.setcon*p.day_to_sec;
    b_ds.b4 = b_ds.b4*p.setcon*p.day_to_sec;
    b_ds.b5 = b_ds.b5*p.setcon*p.day_to_sec;
    b_ds.b25 = b_ds.b25*p.setcon*p.day_to_sec;

    % Pack up the betas
    p2.b1 = b_brown.b1 + b_shear.b1 + b_ds.b1;
    p2.b2 = b_brown.b2 + b_shear.b2 + b_ds.b2;
    p2.b3 = b_brown.b3 + b_shear.b3 + b_ds.b3;
    p2.b4 = b_brown.b4 + b_shear.b4 + b_ds.b4;
    p2.b5 = b_brown.b5 + b_shear.b5 + b_ds.b5;
    p2.b25 = p2.b2 - p2.b3 - p2.b4 - p2.b5;

    % Calculate linear terms
    p2.growth = CalcGrowth(p);
    p2.sink_loss = CalcSinkingLoss(p);
    p2.linear = p2.growth - p2.sink_loss;

    % Calculate disaggregation terms
    p2.disagg_minus = p.c3*diag(p.c4.^(1 : p.n_sections));
    p2.disagg_plus = p.c3*diag(p.c4.^(2:p.n_sections),-1);

    % Initial spectrum
    spec_init_legacy = CalcInitialSpec(p, p2);

    % Solve ODEs
    calcomp = 1:p.n_sections;
    abs_tol = 1.0e-18;
    rel_tol = 3.0e-14;
    at = (abs_tol * 1.5 .^ (-(calcomp-1)));
    t_span = p.t_init : p.delta_t : p.t_final - 1;
    ode_options = odeset('RelTol', rel_tol, 'Refine', 0, 'AbsTol', at, 'Jacobian', @CalcCoagJac);
    [t_out_legacy, y_legacy] = ode15s(@CalcCoagDeriv, t_span, spec_init_legacy, ode_options, p2);

    % Generate legacy plots
    fprintf('  Generating legacy plots...\n');
    outflag_legacy = CoagOutput(p, p2, t_out_legacy, y_legacy);

    fprintf('  Legacy version completed successfully.\n');
    legacy_success = true;
catch ME
    fprintf('  Legacy version failed: %s\n', ME.message);
    legacy_success = false;
end

%% Run OOP Version
fprintf('\nRunning OOP Version...\n');
try
    % Create simulation
    sim = CoagulationSimulation();
    result_oop = sim.run();

    % Generate OOP plots
    fprintf('  Generating OOP plots...\n');
    sim.generateOutputs(true);

    fprintf('  OOP version completed successfully.\n');
    oop_success = true;
catch ME
    fprintf('  OOP version failed: %s\n', ME.message);
    oop_success = false;
end

%% Compare Results
if legacy_success && oop_success
    fprintf('\n=== Detailed Comparison Results ===\n');

    % Compare concentration matrices
    conc_diff = abs(y_legacy - result_oop.concentrations);
    max_conc_diff = max(conc_diff(:));
    mean_conc_diff = mean(conc_diff(:));
    fprintf('Concentration matrices:\n');
    fprintf('  Max difference: %.2e\n', max_conc_diff);
    fprintf('  Mean difference: %.2e\n', mean_conc_diff);
    fprintf('  Relative max diff: %.2e\n', max_conc_diff / max(y_legacy(:)));

    % Check for negative values
    legacy_neg = sum(y_legacy(:) < 0);
    oop_neg = sum(result_oop.concentrations(:) < 0);
    fprintf('\nNegative values:\n');
    fprintf('  Legacy: %d negative values\n', legacy_neg);
    fprintf('  OOP:    %d negative values\n', oop_neg);

    if legacy_neg > 0
        fprintf('  Legacy min: %.2e\n', min(y_legacy(:)));
        fprintf('  Legacy max negative: %.2e\n', max(y_legacy(y_legacy < 0)));
    end
    if oop_neg > 0
        fprintf('  OOP min: %.2e\n', min(result_oop.concentrations(:)));
        fprintf('  OOP max negative: %.2e\n', max(result_oop.concentrations(result_oop.concentrations < 0)));
    end

    % Compare specific time points
    fprintf('\nTime point comparison:\n');
    for i = [1, 5, 10, 15, 20, 25, 30]
        if i <= size(y_legacy, 1) && i <= size(result_oop.concentrations, 1)
            diff_i = abs(y_legacy(i, :) - result_oop.concentrations(i, :));
            fprintf('  Time point %d: max diff = %.2e\n', i, max(diff_i));
        end
    end

    % Save comparison data
    comparison_data = struct();
    comparison_data.legacy_time = t_out_legacy;
    comparison_data.legacy_concentrations = y_legacy;
    comparison_data.oop_time = result_oop.time;
    comparison_data.oop_concentrations = result_oop.concentrations;
    comparison_data.concentration_diff = conc_diff;
    save('detailed_comparison_data.mat', 'comparison_data');
    fprintf('\nComparison data saved to detailed_comparison_data.mat\n');

else
    fprintf('\nCannot compare - one or both versions failed.\n');
end

fprintf('\n=== Plot Comparison ===\n');
fprintf('Check the generated figures to compare:\n');
fprintf('- Figure 1: Spectra and fluxes\n');
fprintf('- Figure 2: Mass balance\n');
fprintf('- Figure 3: Coagulation vs settling losses ratio\n');
fprintf('- Figure 4: 3D surface plots\n');
