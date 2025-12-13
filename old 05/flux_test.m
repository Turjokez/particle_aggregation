function flux_test
% FLUX_TEST
% Adrian-style mass / flux tests on the 1-D column model.
%
% TEST 1: NPP ON, constant conditions → steady state
%   • Run with 3 different "mu" scaling values for NPP
%   • Compare bottom mass flux vs NPP input
%
% TEST 2: NPP OFF, no biology
%   • Initial condition = particles only in surface layer
%   • Check mass budget: M0 ≈ Mf + Exported

    clearvars -except ans; %#ok<CLVAR>
    close all;
    clc;

    fprintf('=====================================================\n');
    fprintf('   ADRIAN MASS / FLUX TESTS (1-D COLUMN)\n');
    fprintf('=====================================================\n\n');

    % --------------------------------------------------------------
    % Common settings
    % --------------------------------------------------------------
    Ns           = 35;        % number of size bins
    z_max        = 200;       % [m]
    dz           = 10;        % [m]
    epsilon_val  = 1e-7;      % [m^2 s^-3] turbulence
    base_NPP     = 1e-7;      % base NPP rate [#/cm^3/day], scaled by mu

    t_final  = 60;            % [days]
    tspan    = [0 t_final];   % IMPORTANT: only [t0 tf] to avoid restart bug

    mu_list  = [0.10 0.50 1.00];

    F_ss          = zeros(size(mu_list));  % steady-state bottom mass flux
    PP_ss         = zeros(size(mu_list));  % steady-state NPP mass input rate
    ratio_F_to_PP = zeros(size(mu_list));  % flux / input

    % ==============================================================
    % TEST 1: NPP STEADY STATE
    % ==============================================================
    fprintf('==================== TEST 1: NPP STEADY STATE ====================\n\n');

    for im = 1:numel(mu_list)
        mu_val = mu_list(im);
        fprintf('--- TEST 1.%d with mu = %.3f d^-1 ---\n', im, mu_val);

        % ------------- Build config for this run -------------------
        cfg = SimulationConfig();
        cfg.use_column           = true;
        cfg.n_sections           = Ns;
        cfg.z_max                = z_max;
        cfg.dz                   = dz;
        cfg.epsilon              = epsilon_val;
        cfg.use_NPP              = true;
        cfg.NPP_rate             = base_NPP * mu_val;
        cfg.attenuation_rate     = 0.0;
        cfg.loss_rate            = 0.0;
        cfg.disagg_use_nonlinear = false;

        % ------------- Run full column model ----------------------
        sim = CoagulationSimulation(cfg);

        % (Optional) relax solver tolerances a bit
        if ~isempty(sim.solver.options)
            sim.solver.options = odeset(sim.solver.options, ...
                'RelTol', 1e-3, 'AbsTol', 1e-8);
        end

        result  = sim.run('tspan', tspan);
        gridObj = sim.grid;      % NOTE: NOT called "grid"

        t  = result.time(:);               % [Nt x 1]
        Y  = result.concentrations;        % [Nt x (Ns*Nz)]
        Nt = size(Y,1);
        Nz = floor(cfg.z_max / cfg.dz);

        % Reshape Y → [Ns x Nz x Nt]
        Y_resh = reshape(Y.', cfg.n_sections, Nz, Nt);

        % Size-bin properties
        av_vol = gridObj.av_vol(:);           % [Ns x 1], volume per particle [cm^3]
        dz_cm  = cfg.dz * 100;                % [cm]

        % Settling velocities [cm s^-1]
        r_cm = gridObj.getFractalRadii();
        r_v  = gridObj.getConservedRadii();
        ws   = SettlingVelocityService.velocity(r_cm, r_v, gridObj.setcon);

        % ---------- Bottom mass flux (per area, per second) -------
        Y_bottom = squeeze(Y_resh(:, Nz, :));       % [Ns x Nt]
        F_mass   = zeros(Nt,1);                    % [cm^3/cm^2/s]

        for k = 1:Nt
            F_mass(k) = sum( ws .* (Y_bottom(:,k) .* av_vol) );
        end

        % Approx steady state → last 20% of available times
        idx_ss   = t > 0.8 * t(end);
        F_ss(im) = mean(F_mass(idx_ss));

        % ---------- NPP mass input (per area, per second) ----------
        % Only smallest bin, top layer
        dz_cm_top = dz_cm;
        G_day     = cfg.NPP_rate * av_vol(1) * dz_cm_top;   % [cm^3/cm^2/day]
        PP_ss(im) = G_day / 86400;                         % [cm^3/cm^2/s]

        ratio_F_to_PP(im) = F_ss(im) / PP_ss(im);

        fprintf('  <F_bottom>_SS (mass flux, s^-1)  = %.3e\n', F_ss(im));
        fprintf('  <PP input>_SS (mass src,   s^-1) = %.3e\n', PP_ss(im));
        fprintf('  Ratio F_out / PP = %.3f\n\n', ratio_F_to_PP(im));
    end

    % ---- Simple diagnostic plot for Test 1 ------------------------
    figure;
    plot(mu_list, ratio_F_to_PP, 'o-', 'LineWidth', 1.5, 'MarkerSize', 8);
    xlabel('\mu (d^{-1})');
    ylabel('F_{bottom} / PP');
    title('Test 1: Steady-state bottom flux vs NPP input');
    grid on;

    % ==============================================================
    % TEST 2: NPP OFF, MASS BUDGET
    % ==============================================================
    fprintf('==================== TEST 2: NPP OFF MASS BUDGET ====================\n\n');

    cfg2 = SimulationConfig();
    cfg2.use_column           = true;
    cfg2.n_sections           = Ns;
    cfg2.z_max                = z_max;
    cfg2.dz                   = dz;
    cfg2.epsilon              = epsilon_val;
    cfg2.use_NPP              = false;
    cfg2.attenuation_rate     = 0.0;
    cfg2.loss_rate            = 0.0;
    cfg2.disagg_use_nonlinear = false;

    gridObj2 = cfg2.derive();   % derived grid for initial spectrum

    % ---------- Initial condition: particles only in surface layer ----------
    N0_slab = InitialSpectrumBuilder.initialSpectrum(cfg2, gridObj2);  % [Ns x 1]

    Nz2      = floor(cfg2.z_max / cfg2.dz);
    N0_grid  = zeros(cfg2.n_sections, Nz2);   % [Ns x Nz]
    N0_grid(:,1) = N0_slab;                   % all particles in top layer

    v0_col = N0_grid.';    % [Nz x Ns]
    v0_col = v0_col(:);    % [Nz*Ns x 1]

    dz_cm2  = cfg2.dz * 100;
    av_vol2 = gridObj2.av_vol(:);

    % Total initial mass per area (biovolume)
    mass0 = dz_cm2 * sum( N0_slab(:) .* av_vol2(:) );
    fprintf('Initial mass M0 (per area) = %.6e [cm^3/cm^2]\n', mass0);

    % ---------- Run the 1-D column with NPP off -------------------
    sim2 = CoagulationSimulation(cfg2);

    if ~isempty(sim2.solver.options)
        sim2.solver.options = odeset(sim2.solver.options, ...
            'RelTol', 1e-3, 'AbsTol', 1e-8);
    end

    result2 = sim2.run('tspan', tspan, 'v0', v0_col);

    t2  = result2.time(:);                  % [Nt2 x 1]
    Y2  = result2.concentrations;           % [Nt2 x (Ns*Nz2)]
    Nt2 = size(Y2,1);
    Y2_resh = reshape(Y2.', cfg2.n_sections, Nz2, Nt2);   % [Ns x Nz x Nt2]

    % Column-integrated mass vs time: ∑_z ∑_j N(z,j)*V_j * dz
    mass_t = zeros(Nt2,1);
    for k = 1:Nt2
        slice     = Y2_resh(:,:,k);                % [Ns x Nz]
        slice_vec = slice(:);                      % [Ns*Nz x 1]
        av_rep    = repmat(av_vol2, Nz2, 1);       % [Ns*Nz x 1]
        mass_t(k) = dz_cm2 * sum( slice_vec .* av_rep );
    end

    % Bottom mass flux vs time
    r_cm2 = gridObj2.getFractalRadii();
    r_v2  = gridObj2.getConservedRadii();
    ws2   = SettlingVelocityService.velocity(r_cm2, r_v2, gridObj2.setcon);

    Y2_bottom = squeeze(Y2_resh(:, Nz2, :));      % [Ns x Nt2]
    F2_mass   = zeros(Nt2,1);
    for k = 1:Nt2
        F2_mass(k) = sum( ws2 .* (Y2_bottom(:,k) .* av_vol2) );
    end

    % Time-integrated exported mass: ∫ F2_mass dt (trapezoidal)
    dt_sec = diff(t2) * 86400;                    % [s]
    F_mid  = 0.5 * (F2_mass(1:end-1) + F2_mass(2:end));
    Export_mass = sum(F_mid .* dt_sec);           % [cm^3/cm^2]

    mass_end = mass_t(end);

    fprintf('Final   mass Mf (per area) = %.6e [cm^3/cm^2]\n', mass_end);
    fprintf('Exported mass E            = %.6e [cm^3/cm^2]\n', Export_mass);
    fprintf('Mf + E - M0                = %.6e (should be ~0)\n', ...
            mass_end + Export_mass - mass0);
    rel_err = (mass_end + Export_mass - mass0) / max(mass0, eps);
    fprintf('Relative error             = %.3e\n', rel_err);

    % ---------- Quick plots for Test 2 ----------------------------
    figure;
    subplot(2,1,1);
    plot(t2, mass_t, 'LineWidth', 1.5);
    xlabel('Time (days)');
    ylabel('Column mass (cm^3/cm^2)');
    title('Test 2: Column-integrated mass (NPP off)');
    grid on;

    subplot(2,1,2);
    plot(t2, F2_mass, 'LineWidth', 1.5);
    xlabel('Time (days)');
    ylabel('Bottom mass flux (cm^3/cm^2/s)');
    title('Test 2: Bottom mass flux');
    grid on;

end
