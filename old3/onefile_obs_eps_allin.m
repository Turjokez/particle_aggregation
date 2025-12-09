function onefile_obs_eps_allin
    fprintf('=== onefile_obs_eps_allin [CORRECTED] ===\n');

    % --- 1. Load Data ---
    matfile = fullfile(pwd,'epsilon_daily.mat');
    if ~isfile(matfile)
        error('epsilon_daily.mat not found in %s', pwd);
    end

    S = load(matfile);
    [t_days, eps_raw, t_key, e_key] = grab_time_eps(S);

    fprintf('MAT file: %s\n', matfile);
    fprintf('Range   : [%.2e, %.2e] W/kg\n', min(eps_raw), max(eps_raw));

    % --- 2. Setup Output Directory ---
    outdir = fullfile(getenv('HOME'), 'Desktop', 'run_obs_eps_nonlinear', ...
        datestr(now,'yyyymmdd_HHMMSS'));
    if ~exist(outdir,'dir'), mkdir(outdir); end
    
    % --- 3. Run Simulations (Coarse vs Fine) ---
    T_total = t_days(end) - t_days(1);
    
    % Run 1: Coarse (Fast)
    N1 = 50;
    t1 = linspace(0, T_total, N1).';
    eps1 = interp1(t_days, eps_raw, t1, 'linear', 'extrap');
    
    fprintf('\n--- Running SIM 1 (N=%d) ---\n', N1);
    sim1 = run_sim_with_eps(t1, eps1);
    out1 = sim1.exportMinimalOutputs();
    
    % === FIX FOR BLANK FIGURES ===
    % Manually calculate diameter if the class didn't grab it
    if isempty(out1.D_um)
        % Get radii in cm, convert to diameter in microns
        r_cm = sim1.grid.getFractalRadii(); 
        out1.D_um = 2 * r_cm(:) * 1e4; 
    end

    % Run 2: Fine (Better Physics)
    N2 = 100;
    t2 = linspace(0, T_total, N2).';
    eps2 = interp1(t_days, eps_raw, t2, 'linear', 'extrap');
    
    fprintf('\n--- Running SIM 2 (N=%d) ---\n', N2);
    sim2 = run_sim_with_eps(t2, eps2);
    out2 = sim2.exportMinimalOutputs();
    
    % === FIX FOR BLANK FIGURES ===
    if isempty(out2.D_um)
        r_cm = sim2.grid.getFractalRadii(); 
        out2.D_um = 2 * r_cm(:) * 1e4; 
    end

    % --- 4. Generate Plots ---
    fprintf('\nGenerating plots in %s...\n', outdir);
    make_eps_export_plot(t1, eps1, out1, fullfile(outdir,'observed_eps_nonlinear.png'));
    make_dt_robustness_plot(out1, out2, fullfile(outdir,'dt_robustness.png'));
    make_mass_vs_settling_plot(out1, fullfile(outdir,'mass_vs_settling.png'));

    fprintf('✓ All done.\n');
end

% =====================================================================
% CONFIGURATION & RUNNER
% =====================================================================
function sim = run_sim_with_eps(t_days, eps_series)
    cfg         = SimulationConfig();
    cfg.t_init  = t_days(1);
    cfg.t_final = t_days(end);
    cfg.delta_t = mean(diff(t_days));
    
    % === PHYSICS SWITCHES ===
    cfg.disagg_use_nonlinear = true; 
    cfg.disagg_kmax_a         = 0.70;  
    cfg.disagg_beta           = 0.35;  
    
    % Primary Production (Keeps mass in system)
    cfg.use_NPP     = true;
    cfg.NPP_rate    = 2e-4;       
    cfg.NPP_profile = 'constant'; 

    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_days, eps_series);
    
    % Run
    sim.run('tspan', t_days); 
end

% =====================================================================
% DATA LOADER HELPER (FIXED TYPO)
% =====================================================================
function [t_days, eps_vec, time_key, eps_key] = grab_time_eps(S)

    fn = fieldnames(S);
    while numel(fn) == 1 && isstruct(S.(fn{1}))
        S  = S.(fn{1});
        fn = fieldnames(S);
    end

    numNames = {};
    numVals  = {};
    for k = 1:numel(fn)
        v = S.(fn{k});
        if isnumeric(v) && numel(v) > 1
            numNames{end+1} = fn{k};
            numVals{end+1}  = v;
        end
    end

    if isempty(numNames)
        error('grab_time_eps: MAT-file has no numeric fields.');
    end

    eps_key = '';
    for k = 1:numel(numNames)
        lname = lower(numNames{k});
        % === TYPO FIXED HERE (lame -> lname) ===
        if contains(lname,'eps') || contains(lname,'epsilon') || contains(lname,'dissip')
            eps_key = numNames{k};
            break
        end
    end
    if isempty(eps_key)
        lens = cellfun(@numel,numVals);
        [~,ik] = max(lens);
        eps_key = numNames{ik};
    end

    eps_vec = S.(eps_key);
    eps_vec = eps_vec(:);

    N = numel(eps_vec);
    time_key = '';
    for k = 1:numel(fn)
        v = S.(fn{k});
        if ~isnumeric(v) || numel(v) ~= N, continue; end
        lname = lower(fn{k});
        if contains(lname,'time') || contains(lname,'mt') || ...
           contains(lname,'day')  || contains(lname,'date') || strcmpi(lname,'t')
            time_key = fn{k};
            break
        end
    end
    if isempty(time_key)
        for k = 1:numel(numNames)
            if numel(numVals{k}) == N && ~strcmp(numNames{k},eps_key)
                time_key = numNames{k};
                break
            end
        end
    end

    if isempty(time_key)
        time_key = 'synthetic';
        t_days   = (0:N-1).';
    else
        t_raw = S.(time_key);
        t_raw = t_raw(:);
        if max(t_raw) > 1e4 && range(t_raw) > 0
            t_days = t_raw - t_raw(1);
        else
            t_days = t_raw;
        end
    end

    mask   = isfinite(t_days) & isfinite(eps_vec);
    t_days = t_days(mask);
    eps_vec = eps_vec(mask);

    [t_days, idx] = sort(t_days(:));
    eps_vec = eps_vec(idx);
    [t_days, ia] = unique(t_days, 'stable');
    eps_vec = eps_vec(ia);

    t0     = t_days(1);
    t_days = t_days - t0;
    eps_vec = max(eps_vec, 1e-12);
end

% =====================================================================
% PLOTTING HELPERS
% =====================================================================
function make_eps_export_plot(t_days, eps_series, out, fname)
    f = figure('Color','w','Position',[100 100 1200 750]);

    % Top Left: Epsilon Input
    subplot(2,2,1);
    semilogy(t_days, eps_series, 'k-','LineWidth',1.2);
    yline(1e-6, 'r--', 'Threshold');
    grid on
    xlabel('Time (days)');
    ylabel('\epsilon (W kg^{-1})');
    title('Input Turbulence');

    % Top Right: Export Flux Proxy
    subplot(2,2,2);
    t  = out.t_days(:);
    flux_inst = out.settling_loss(:); 
    plot(t, flux_inst, 'b-','LineWidth',1.2);
    grid on
    xlabel('Time (days)');
    ylabel('Instant Flux (Model Units)');
    title('Particle Export Flux');

    % Bottom Left: PSD comparison (Start vs End)
    subplot(2,2,3);
    if ~isempty(out.N) && ~isempty(out.D_um)
        loglog(out.D_um(:), max(out.N(:,1),1e-20),'b-','LineWidth',1.0); hold on
        loglog(out.D_um(:), max(out.N(:,end),1e-20),'r-','LineWidth',1.0); hold off
        grid on
        legend('Start','End');
        xlabel('Diameter (\mum)');
        ylabel('# cm^{-3}');
        title('PSD Snapshot');
    else
        text(0.5,0.5,'No PSD Data');
    end

    % Bottom Right: The "Flip" Heatmap
    subplot(2,2,4);
    if ~isempty(out.N) && ~isempty(out.D_um)
        imagesc(t, log10(out.D_um(:)), log10(max(out.N,1e-20)));
        set(gca,'YDir','normal');
        cb = colorbar;
        cb.Label.String = 'log_{10} Concentration';
        colormap(jet); 
        caxis([-5 0]); 
        xlabel('Time (days)');
        ylabel('log_{10} D (\mum)');
        title('PSD Evolution (Look for Dip during Storms)');
    else
        text(0.5,0.5,'No PSD Data');
    end

    exportgraphics(f,fname,'Resolution',150);
    close(f);
end

function make_dt_robustness_plot(out1, out2, fname)
    t1  = out1.t_days(:);
    t2  = out2.t_days(:);
    cum1 = cumsum(out1.settling_loss(:)); 
    cum2 = cumsum(out2.settling_loss(:)); 
    cum1 = cum1 / max(cum1);
    cum2 = cum2 / max(cum2);

    f = figure('Color','w','Position',[50 50 1000 400]);
    plot(t1, cum1, 'b-','LineWidth',1.5); hold on
    plot(t2, cum2, 'r--','LineWidth',1.5);
    grid on
    xlabel('Time (days)');
    ylabel('Normalized Cumulative Flux');
    title('Time-Step Robustness Check');
    legend('Coarse Grid','Fine Grid','Location','SouthEast');
    exportgraphics(f,fname,'Resolution',150);
    close(f);
end

function make_mass_vs_settling_plot(out, fname)
    t  = out.t_days(:);
    rate_sett = out.settling_loss(:); 
    f = figure('Color','w','Position',[50 50 1000 400]);
    yyaxis left
    plot(t, rate_sett, 'b-','LineWidth',1.5);
    ylabel('Settling Loss Rate');
    yyaxis right
    mass_total = sum(out.N, 1);
    plot(t, mass_total, 'k-','LineWidth',1.0);
    ylabel('Total Suspended Mass');
    grid on
    xlabel('Time (days)');
    title('Dynamics Check: Settling vs Total Mass');
    exportgraphics(f,fname,'Resolution',150);
    close(f);
end