%% Final Verification: Legacy vs OOP Plot Comparison
% This script generates both versions and creates side-by-side plots for comparison

clear all; close all;

fprintf('=== Final Verification: Legacy vs OOP Plot Comparison ===\n\n');

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

    % Generate legacy output data
    legacy_data = generateLegacyOutputData(p, p2, t_out_legacy, y_legacy);

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

    % Generate OOP output data
    oop_data = generateOOPOutputData(result_oop, sim.grid, sim.config);

    fprintf('  OOP version completed successfully.\n');
    oop_success = true;
catch ME
    fprintf('  OOP version failed: %s\n', ME.message);
    oop_success = false;
end

%% Generate Comparison Plots
if legacy_success && oop_success
    fprintf('\nGenerating comparison plots...\n');

    % Figure 1: Number spectrum comparison
    figure(1);
    clf;

    % Legacy plot (left)
    subplot(2, 2, 1);
    r_i = p.amfrac * p.av_vol.^p.bmfrac;
    r_v = (0.75/pi*p.av_vol).^(1.0/3.0);
    diam_i = 2.0*p.r_to_rg*r_i;
    diaratio = (p.fr_dim/3) * r_v ./ r_i;
    nspec_init = y_legacy(1, :) ./ (1.5*p.v_lower') ./ p.dwidth';
    nspec_final = y_legacy(end, :) ./ (1.5*p.v_lower') ./ p.dwidth';
    ns_init = nspec_init .* diaratio;
    ns_final = nspec_final .* diaratio;
    loglog(diam_i, ns_init, 'b', diam_i, ns_final, 'r');
    xlabel('Particle diameter [cm]');
    ylabel('Number spectrum [# cm^{-4}]');
    title('Legacy: Number Spectrum');
    axis tight;

    % OOP plot (right)
    subplot(2, 2, 2);
    loglog(oop_data.diam_i, oop_data.nspec_i(1, :), 'b', ...
        oop_data.diam_i, oop_data.nspec_i(end, :), 'r');
    xlabel('Particle diameter [cm]');
    ylabel('Number spectrum [# cm^{-4}]');
    title('OOP: Number Spectrum');
    axis tight;

    % Legacy sectional concentration
    subplot(2, 2, 3);
    semilogy(t_out_legacy, y_legacy, t_out_legacy, sum(y_legacy,2), '*--');
    xlabel('Time [d]');
    ylabel('Sectional concentration [vol/vol/sect]');
    title('Legacy: Sectional Concentration');

    % OOP sectional concentration
    subplot(2, 2, 4);
    semilogy(result_oop.time, result_oop.concentrations, ...
        result_oop.time, sum(result_oop.concentrations,2), '*--');
    xlabel('Time [d]');
    ylabel('Sectional concentration [vol/vol/sect]');
    title('OOP: Sectional Concentration');

    % Figure 2: Coagulation vs Settling Losses Ratio
    figure(2);
    clf;

    % Calculate ratios for both versions
    try
        [sec_gains_legacy, sec_losses_legacy] = SectionalMassBalance(y_legacy, p2);
        [total_gains_legacy, total_losses_legacy] = TotalMassBalance(y_legacy, p2);
    catch
        % Fallback calculation
        total_losses_legacy = struct();
        total_losses_legacy.coag = zeros(size(t_out_legacy));
        total_losses_legacy.sett = zeros(size(t_out_legacy));
    end

    subplot(1, 2, 1);
    ratio_legacy = total_losses_legacy.coag ./ total_losses_legacy.sett;
    plot(t_out_legacy, ratio_legacy);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 14);
    xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 14);
    ylabel('(Coag Losses)/(Settling Losses)', 'FontName', 'Helvetica', 'FontSize', 14);
    title('Legacy: Coag vs Sett Ratio');

    subplot(1, 2, 2);
    ratio_oop = result_oop.diagnostics.total_losses.coag ./ result_oop.diagnostics.total_losses.sett;
    plot(result_oop.time, ratio_oop);
    set(gca, 'FontName', 'Helvetica', 'FontSize', 14);
    xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 14);
    ylabel('(Coag Losses)/(Settling Losses)', 'FontName', 'Helvetica', 'FontSize', 14);
    title('OOP: Coag vs Sett Ratio');

    % Figure 3: Direct comparison of key metrics
    figure(3);
    clf;

    % Total mass comparison
    subplot(2, 2, 1);
    plot(t_out_legacy, sum(y_legacy,2), 'b-', 'LineWidth', 2);
    hold on;
    plot(result_oop.time, sum(result_oop.concentrations,2), 'r--', 'LineWidth', 2);
    xlabel('Time [d]');
    ylabel('Total Mass');
    title('Total Mass Comparison');
    legend('Legacy', 'OOP', 'Location', 'best');
    grid on;

    % Difference in total mass
    subplot(2, 2, 2);
    mass_diff = abs(sum(y_legacy,2) - sum(result_oop.concentrations,2));
    semilogy(t_out_legacy, mass_diff);
    xlabel('Time [d]');
    ylabel('|Mass Difference|');
    title('Total Mass Difference');
    grid on;

    % Coagulation ratio comparison
    subplot(2, 2, 3);
    plot(t_out_legacy, ratio_legacy, 'b-', 'LineWidth', 2);
    hold on;
    plot(result_oop.time, ratio_oop, 'r--', 'LineWidth', 2);
    xlabel('Time [d]');
    ylabel('(Coag Losses)/(Settling Losses)');
    title('Coagulation Ratio Comparison');
    legend('Legacy', 'OOP', 'Location', 'best');
    grid on;

    % Difference in coagulation ratio
    subplot(2, 2, 4);
    ratio_diff = abs(ratio_legacy - ratio_oop);
    semilogy(t_out_legacy, ratio_diff);
    xlabel('Time [d]');
    ylabel('|Ratio Difference|');
    title('Coagulation Ratio Difference');
    grid on;

    fprintf('Comparison plots generated successfully!\n');
    fprintf('Check Figures 1-3 for visual comparison.\n');

    % Print summary statistics
    fprintf('\n=== Summary Statistics ===\n');
    fprintf('Concentration max difference: %.2e\n', max(abs(y_legacy(:) - result_oop.concentrations(:))));
    fprintf('Total mass max difference: %.2e\n', max(abs(sum(y_legacy,2) - sum(result_oop.concentrations,2))));
    fprintf('Coagulation ratio max difference: %.2e\n', max(abs(ratio_legacy - ratio_oop)));

else
    fprintf('\nCannot generate comparison plots - one or both versions failed.\n');
end

%% Helper Functions
function data = generateLegacyOutputData(p, p2, t_out, y)
% Generate output data similar to CoagOutput.m
data = struct();
data.t = t_out;
data.Y = y;

% Calculate additional spectra
n_times = length(t_out);
nspec_v = zeros(n_times, p.n_sections);
masspec_v = nspec_v;
fluxsect = nspec_v;
fluxspec = nspec_v;

r_i = p.amfrac * p.av_vol.^p.bmfrac;
r_v = (0.75/pi*p.av_vol).^(1.0/3.0);
set_vel = SettlingVelocity(r_i, r_v, p.setcon);
set_vel = set_vel/100*p.day_to_sec;
diam_i = 2.0*p.r_to_rg*r_i;
diam_v = 2.0*r_v;

for jindx = 1:n_times
    yout = y(jindx,:);
    nspec_v(jindx,:) = yout./(1.5*p.v_lower')./p.dwidth';
    masspec_v(jindx,:) = yout./p.dwidth';
    fluxsect(jindx,:) = yout.*set_vel'*1e6;
    fluxspec(jindx,:) = masspec_v(jindx,:).*set_vel'*1e6;
end

data.nspec_v = nspec_v;
data.masspec_v = masspec_v;
data.fluxsect = fluxsect;
data.fluxspec = fluxspec;
data.diam_i = diam_i;
data.diam_v = diam_v;
data.set_vel = set_vel;
data.v_lower = p.v_lower;
data.dwidth = p.dwidth;

% Calculate image-based spectra
diaratio = (p.fr_dim/3) * diam_v ./ diam_i;
data.nspec_i = nspec_v .* diaratio;
data.masspec_i = masspec_v .* diaratio;
data.fluxspec_i = fluxspec .* diaratio;
end

function data = generateOOPOutputData(result, grid, config)
% Generate output data from OOP results
data = struct();
data.t = result.time;
data.Y = result.concentrations;

% Use the precomputed output data if available
if isfield(result, 'output_data')
    data.nspec_v = result.output_data.nspec_v;
    data.masspec_v = result.output_data.masspec_v;
    data.fluxsect = result.output_data.fluxsect;
    data.fluxspec = result.output_data.fluxspec;
    data.diam_i = result.output_data.diam_i;
    data.diam_v = result.output_data.diam_v;
    data.set_vel = result.output_data.set_vel;
    data.v_lower = result.output_data.v_lower;
    data.dwidth = result.output_data.dwidth;
    data.nspec_i = result.output_data.nspec_i;
    data.masspec_i = result.output_data.masspec_i;
    data.fluxspec_i = result.output_data.fluxspec_i;
else
    % Calculate manually if not available
    data = generateLegacyOutputData(config, struct(), result.time, result.concentrations);
end
end
