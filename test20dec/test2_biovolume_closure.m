function test2_biovolume_closure()
% TEST 2: Column closure diagnostic (BIOVOLUME/AREA)
% v treated as biovolume concentration (cm^3/cm^3)

    fprintf('TEST 2: Column closure diagnostic (BIO-VOLUME)\n');

    run_case('A (coag OFF, disagg OFF)', false);
    run_case('B (coag ON,  disagg OFF)', true);
end

function run_case(label, coag_on)

    fprintf('TEST 2%s: Column closure diagnostic (BIO-VOLUME)\n', label);

    cfg = SimulationConfig();
    cfg.use_column = true;

    cfg.z_max = 65;
    cfg.dz    = 5;

    % PURE SINKING TEST:
    cfg.growth = 0;
    cfg.growth_mode = 'shift';

    % *** KEY FIX: TURN OFF DISAGGREGATION ***
    if isprop(cfg,'c3'); cfg.c3 = 0; end

    % Flux-form sinking
    cfg.sinking_form = 'flux';

    % Coag toggle
    if isprop(cfg,'enable_coag')
        cfg.enable_coag = logical(coag_on);
    end

    % time
    if isprop(cfg,'t_init');  cfg.t_init  = 0;   end
    if isprop(cfg,'t_final'); cfg.t_final = 120; end
    if isprop(cfg,'delta_t'); cfg.delta_t = 1;   end

    sim = CoagulationSimulation(cfg);
    out = sim.run();

    t = out.time(:);
    Y = out.concentrations;

    Ns = cfg.n_sections;
    Nz = cfg.getNumLayers();
    dz_cm = cfg.dz * 100;

    [S, ~, w_cmday] = LinearProcessBuilder.sinkingMatrix(cfg, sim.grid);
    w_cmday = w_cmday(:);

    Nt = numel(t);

    Bcol    = zeros(Nt,1);  % cm^3/cm^2
    Fbottom = zeros(Nt,1);  % cm^3/cm^2/day
    dBdt_op = zeros(Nt,1);  % cm^3/cm^2/day (from sinking operator)

    for it = 1:Nt
        v = Y(it,:)';
        V = reshape(v,[Ns Nz]);

        Bcol(it) = dz_cm * sum(v);

        vbot = V(:,end);
        Fbottom(it) = sum(w_cmday .* vbot);

        dv_sink = -(S*v);
        dBdt_op(it) = dz_cm * sum(dv_sink);
    end

    dBdt_ts = gradient(Bcol, t);
    intF    = cumtrapz(t, Fbottom);

    closure = Bcol + intF;
    closure_err_end = closure(end) - Bcol(1);

    RelErr_operator_rate = max(abs(dBdt_op + Fbottom)) / max(max(abs(Fbottom)), eps);
    RelErr_timeseries_vs_operator = max(abs(dBdt_ts - dBdt_op)) / max(max(abs(dBdt_ts)), eps);
    RelErr_integral = abs(closure_err_end) / max(abs(Bcol(1)), eps);

    fprintf('\n[BIOVOLUME/AREA UNITS]\n');
    fprintf('Initial B0                   = %.6e  [cm^3/cm^2]\n', Bcol(1));
    fprintf('Final   B(end)               = %.6e  [cm^3/cm^2]\n', Bcol(end));
    fprintf('Final   Int Export dt         = %.6e  [cm^3/cm^2]\n', intF(end));
    fprintf('Final   closure(end)          = %.6e  [cm^3/cm^2]\n', closure(end));
    fprintf('Final   closure error(end)    = %.6e  [cm^3/cm^2]\n', closure_err_end);

    fprintf('\nRelErr(operator rate)         = %.6e\n', RelErr_operator_rate);
    fprintf('RelErr(timeseries vs operator)= %.6e\n', RelErr_timeseries_vs_operator);
    fprintf('RelErr(integral)              = %.6e\n', RelErr_integral);

    % ================= FIGURES =================
    fig1 = figure('Name',['Test2 ',label],'Color','w');
    plot(t, Bcol, 'LineWidth', 2); hold on;
    plot(t, closure, '--', 'LineWidth', 2);
    xlabel('time [day]'); ylabel('Bcol, closure [cm^3/cm^2]');
    legend('Bcol(t)', 'Bcol + \int F dt', 'Location','best');
    grid on;
    title(['Closure check: ', label], 'Interpreter','none');

    fig2 = figure('Name',['Test2 rates ',label],'Color','w');
    plot(t, dBdt_ts, 'LineWidth', 2); hold on;
    plot(t, dBdt_op, '--', 'LineWidth', 2);
    plot(t, -Fbottom, ':', 'LineWidth', 2);
    xlabel('time [day]'); ylabel('rate [cm^3/cm^2/day]');
    legend('dB/dt (timeseries)', 'dB/dt (operator sink)', '-Fbottom', 'Location','best');
    grid on;
    title(['Rate consistency: ', label], 'Interpreter','none');

    % Save figures
    safe = matlab.lang.makeValidName(strrep(label,' ','_'));
    saveas(fig1, ['TEST2_closure_', safe, '.png']);
    saveas(fig2, ['TEST2_rates_',   safe, '.png']);

end