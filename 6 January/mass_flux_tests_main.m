function mass_flux_tests_main
% MASS_FLUX_TESTS_MAIN
% Adrian-style mass / flux tests using the CLEAN 0-D coagulation model.
%
% Model here is a single slab (no vertical grid):
%   • "Bottom flux" = flux leaving the slab by settling
%   • "Production"  = growth operator acting on all sections
%
% TEST 1:
%   - Keep everything constant
%   - Run to steady state
%   - Use three growth values (mu = 0.1, 0.5, 1.0 d^-1)
%   - Compare flux leaving the slab to the total mass input from growth
%     (computed with the actual growth matrix)
%
% TEST 2:
%   - Turn growth OFF (no primary production)
%   - Initial condition: all mass placed in the first size class
%   - Run forward
%   - Check mass budget:  M0 ≈ Mf + Exported
%
% IMPORTANT: For these conservation tests we TURN OFF the legacy
%            disaggregation term by setting cfg.c3 = 0.

    clearvars -except ans; %#ok<CLVAR>
    close all;
    clc;

    fprintf('=====================================================\n');
    fprintf('   MASS / FLUX TESTS ON CLEAN MAIN MODEL (0-D SLAB)\n');
    fprintf('=====================================================\n\n');

    % --------------------------------------------------------------
    % COMMON SETTINGS
    % --------------------------------------------------------------
    t_final_days = 60;                 % runtime for both tests [days]
    mu_list      = [0.1, 0.5, 1.0];    % d^-1 growth values (Adrian's mu)

    F_ss   = zeros(size(mu_list));  % steady-state bottom mass flux
    PP_ss  = zeros(size(mu_list));  % steady-state mass input from growth
    ratio  = zeros(size(mu_list));  % F_out / input

    % ==============================================================
    % TEST 1: CONSTANT FORCING, STEADY STATE
    % ==============================================================
    fprintf('==================== TEST 1: STEADY STATE ====================\n\n');

    for im = 1:numel(mu_list)
        mu_val = mu_list(im);
        fprintf('--- TEST 1.%d with mu = %.3f d^-1 ---\n', im, mu_val);

        % ----- Build config for this mu -----
        cfg = SimulationConfig();
        cfg.t_init  = 0.0;
        cfg.t_final = t_final_days;
        cfg.delta_t = 1.0;        % output once per day
        cfg.growth  = mu_val;     % use growth as Adrian's mu

        % Turn OFF legacy disaggregation for conservation test
        cfg.c3 = 0;               % no breakup term

        % ----- Run simulation -----
        sim    = CoagulationSimulation(cfg);
        result = sim.run();       % default initial spectrum

        gridObj = sim.grid;

        % Time and concentrations
        t_days = result.time(:);          % [Nt x 1], days
        Y      = result.concentrations;   % [Nt x Ns]
        Nt = size(Y,1);
        Ns = size(Y,2);

        % Particle and settling properties
        av_vol   = gridObj.av_vol(:);                 % [Ns x 1], cm^3/particle
        dz_cm    = cfg.dz * 100;                      % layer thickness [cm]
        r_cm     = gridObj.getFractalRadii();         % [cm]
        r_v      = gridObj.getConservedRadii();       % [cm]
        ws_cm_s  = SettlingVelocityService.velocity(r_cm, r_v, gridObj.setcon);
        day2sec  = cfg.day_to_sec;                    % 86400

        % Growth matrix (per day) – this is how the model really applies growth
        Gmat = LinearProcessBuilder.growthMatrix(cfg, gridObj);
        % Gmat * N has units [#/cm^3 per day]

        % ---- Compute bottom flux and growth input vs time ----
        F_mass = zeros(Nt,1);   % bottom flux [cm^3 cm^-2 s^-1]
        PP_day = zeros(Nt,1);   % mass input from growth [cm^3 cm^-2 d^-1]

        for k = 1:Nt
            N_k   = Y(k,:).';
            N_pos = max(N_k, 0);   % enforce non-negative concentrations

            % Bottom mass flux (settling out of slab)
            F_mass(k) = sum( ws_cm_s .* (N_pos .* av_vol) );

            % Mass input from growth:
            %   dNdt_growth = Gmat * N_pos  [#/cm^3 / day]
            %   Input_mass  = dz * sum( dNdt_growth .* V )
            dNdt_growth = Gmat * N_pos;
            PP_day(k)   = dz_cm * sum( dNdt_growth .* av_vol );
        end

        % Convert flux to per-day units for comparison
        F_day = F_mass * day2sec;   % [cm^3 cm^-2 d^-1]

        % ---- Steady-state averages (last 20% of record) ----
        t_end  = t_days(end);
        idx_ss = t_days > 0.8 * t_end;
        if ~any(idx_ss)
            idx_ss = true(size(t_days));   % fallback: use all times
        end

        F_ss(im)  = mean(F_day(idx_ss));
        PP_ss(im) = mean(PP_day(idx_ss));
        ratio(im) = F_ss(im) / max(PP_ss(im), eps);

        fprintf('  <F_bottom>_SS (mass flux)  = %.3e [cm^3 cm^-2 d^-1]\n', F_ss(im));
        fprintf('  <Growth input>_SS          = %.3e [cm^3 cm^-2 d^-1]\n', PP_ss(im));
        fprintf('  Ratio F_out / input        = %.3f\n\n', ratio(im));
    end

    % ---- Plot Test 1 diagnostic ----
    figure('Name','Test 1: F_out / Input vs mu','Color','w');
    plot(mu_list, ratio, 'o-','LineWidth',1.8,'MarkerSize',8);
    xlabel('\mu (d^{-1})');
    ylabel('F_{bottom} / Input');
    title('Test 1: Steady-state bottom flux vs growth input');
    grid on;

    % ==============================================================
    % TEST 2: GROWTH OFF, MASS BUDGET
    % ==============================================================
    fprintf('==================== TEST 2: MASS BUDGET ====================\n\n');

    cfg2 = SimulationConfig();
    cfg2.t_init  = 0.0;
    cfg2.t_final = t_final_days;
    cfg2.delta_t = 1.0;
    cfg2.growth  = 0.0;    % growth OFF (no primary production)

    % Turn OFF legacy disaggregation here as well
    cfg2.c3 = 0;

    % Build grid for initial condition and properties
    gridObj2 = cfg2.derive();
    av_vol2  = gridObj2.av_vol(:);
    dz_cm2   = cfg2.dz * 100;
    day2sec2 = cfg2.day_to_sec;

    Ns2 = cfg2.n_sections;

    % ----- Base spectrum (for reference mass) -----
    N0_base = InitialSpectrumBuilder.initialSpectrum(cfg2, gridObj2);  % [Ns2 x 1]

    % Total mass of this base spectrum
    mass_base = dz_cm2 * sum(N0_base .* av_vol2);

    % ----- New initial condition: all mass in first size class -----
    N0 = zeros(Ns2,1);
    if mass_base > 0
        % mass in one particle of bin 1 over the slab:
        mass_per_particle_1 = dz_cm2 * av_vol2(1);
        N0(1) = mass_base / mass_per_particle_1;
    else
        % fallback: just use base N0
        N0 = N0_base;
    end

    % Initial mass per area
    mass0 = dz_cm2 * sum(N0 .* av_vol2);
    fprintf('Initial mass M0 (per area) = %.6e [cm^3 cm^{-2}]\n', mass0);

    % ----- Run simulation with growth off -----
    sim2    = CoagulationSimulation(cfg2);
    result2 = sim2.run('v0', N0);

    t2  = result2.time(:);              % [Nt2 x 1], days
    Y2  = result2.concentrations;       % [Nt2 x Ns2]
    Nt2 = size(Y2,1);

    % ----- Column-integrated mass vs time -----
    mass_t = zeros(Nt2,1);
    for k = 1:Nt2
        N_k = max(Y2(k,:).', 0);        % clamp negatives
        mass_t(k) = dz_cm2 * sum( N_k .* av_vol2 );
    end

    % ----- Bottom mass flux vs time -----
    r_cm2    = gridObj2.getFractalRadii();
    r_v2     = gridObj2.getConservedRadii();
    ws_cm_s2 = SettlingVelocityService.velocity(r_cm2, r_v2, gridObj2.setcon);

    F2_mass = zeros(Nt2,1);             % [cm^3 cm^-2 s^-1]
    for k = 1:Nt2
        N_k = max(Y2(k,:).', 0);
        F2_mass(k) = sum( ws_cm_s2 .* (N_k .* av_vol2) );
    end

    % Time-integrated exported mass (trapezoidal rule)
    if Nt2 > 1
        dt_sec = diff(t2) * day2sec2;            % [s]
        F_mid  = 0.5 * (F2_mass(1:end-1) + F2_mass(2:end));
        Export_mass = sum(F_mid .* dt_sec);      % [cm^3 cm^-2]
    else
        Export_mass = 0;
    end

    mass_end = mass_t(end);

    fprintf('Final   mass Mf (per area) = %.6e [cm^3 cm^{-2}]\n', mass_end);
    fprintf('Exported mass E            = %.6e [cm^3 cm^{-2}]\n', Export_mass);
    fprintf('Mf + E - M0                = %.6e (should be ~0)\n', ...
            mass_end + Export_mass - mass0);

    rel_err = (mass_end + Export_mass - mass0) / max(mass0, eps);
    fprintf('Relative error             = %.3e\n', rel_err);

    % ----- Plots for Test 2 -----
    figure('Name','Test 2: Mass budget','Color','w');

    subplot(2,1,1);
    plot(t2, mass_t, 'LineWidth',1.5);
    xlabel('Time (days)');
    ylabel('Column mass (cm^3 cm^{-2})');
    title('Test 2: Column-integrated mass (growth off, c3=0)');
    grid on;

    subplot(2,1,2);
    plot(t2, F2_mass * day2sec2, 'LineWidth',1.5);
    xlabel('Time (days)');
    ylabel('Bottom mass flux (cm^3 cm^{-2} d^{-1})');
    title('Test 2: Bottom mass flux (growth off, c3=0)');
    grid on;

end
  