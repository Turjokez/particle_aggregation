function test1_steady_state_flux()
% TEST 1 : steady forcing vs bottom export (PP ON)
% v is BIOVOLUME concentration (cm^3/cm^3)

    mu_list = [0.1 0.5 1.0];

    fprintf('TEST 1: Steady-state flux check (PP ON)\n');

    for j = 1:numel(mu_list)
        mu = mu_list(j);

        fprintf('TEST 1.%d  mu = %.3f d^-1\n', j, mu);

        cfg = SimulationConfig();
        cfg.use_column = true;

        cfg.z_max = 65;
        cfg.dz    = 5;
        dz_cm = cfg.dz * 100;

        % PP forcing handled in CoagulationRHS
        cfg.growth_mode = 'pp';
        cfg.growth      = mu;

        % Sinking operator
        cfg.sinking_form = 'flux';

        % Turn OFF disaggregation for this budget test
        if isprop(cfg,'c3'); cfg.c3 = 0; end

        % ---- VALID init_profile + TINY initial condition ----
        if isprop(cfg,'init_profile')
            cfg.init_profile = 'top_only';     % VALID option in your validate()
        end

        tiny = 1e-16;  % try 1e-20 if still too big
        if isprop(cfg,'num_1'); cfg.num_1 = tiny; end
        if isprop(cfg,'n1');    cfg.n1    = tiny; end

        % time
        if isprop(cfg,'t_init');  cfg.t_init  = 0;   end
        if isprop(cfg,'t_final'); cfg.t_final = 400; end
        if isprop(cfg,'delta_t'); cfg.delta_t = 1;   end

        sim = CoagulationSimulation(cfg);
        out = sim.run();

        t = out.time(:);
        Y = out.concentrations;

        Ns = cfg.n_sections;
        Nz = cfg.getNumLayers();

        % w [cm/day]
        [~, ~, w_cmday] = LinearProcessBuilder.sinkingMatrix(cfg, sim.grid);
        w_cmday = w_cmday(:);

        Nt = numel(t);
        Export_flux = zeros(Nt,1);
        Input_flux  = zeros(Nt,1);
        Bcol        = zeros(Nt,1);

        for it = 1:Nt
            v = Y(it,:)';
            V = reshape(v, [Ns Nz]);

            Bcol(it) = dz_cm * sum(v);

            vbot = V(:, end);
            Export_flux(it) = sum(w_cmday .* vbot);

            v1_surface = V(1,1);
            Input_flux(it) = mu * v1_surface * dz_cm;
        end

        % last 30 days mean
        mask = t >= (t(end) - 30);
        mean_export = mean(Export_flux(mask));
        mean_input  = mean(Input_flux(mask));
        ratio_ss = mean_export / max(mean_input, eps);

        % integrated budget
        IntExport = trapz(t, Export_flux);
        IntInput  = trapz(t, Input_flux);
        closure_err = (Bcol(end) + IntExport) - (Bcol(1) + IntInput);

        fprintf('B0 (inventory)               = %.6e  [cm^3/cm^2]\n', Bcol(1));
        fprintf('B(end)                       = %.6e  [cm^3/cm^2]\n', Bcol(end));
        fprintf('Int Export dt                = %.6e  [cm^3/cm^2]\n', IntExport);
        fprintf('Int Input  dt                = %.6e  [cm^3/cm^2]\n', IntInput);
        fprintf('Closure error (end)          = %.6e  [cm^3/cm^2]\n', closure_err);

        fprintf('Mean Export_flux (last ~30d) = %.6e  [cm^3/cm^2/d]\n', mean_export);
        fprintf('Mean Input_PP    (last ~30d) = %.6e  [cm^3/cm^2/d]\n', mean_input);
        fprintf('Ratio Export/Input (SS)      = %.6f\n', ratio_ss);

        % plot
        fig = figure('Color','w');
        semilogy(t, Export_flux, 'LineWidth', 2); hold on;
        semilogy(t, Input_flux,  '--', 'LineWidth', 2);
        xlabel('time [day]'); ylabel('flux [cm^3/cm^2/day]');
        legend('Export bottom','Input PP','Location','best');
        grid on;
        title(sprintf('Test1 mu=%.3f: Export vs Input', mu));
        saveas(fig, sprintf('TEST1_mu_%0.3f_flux.png', mu));
    end
end