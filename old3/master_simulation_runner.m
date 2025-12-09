function master_simulation_runner
    % MASTER_SIMULATION_RUNNER
    % All-in-one script for Particle Aggregation Project.
    
    clc;
    fprintf('======================================================\n');
    fprintf('   PARTICLE AGGREGATION MODEL - MASTER RUNNER\n');
    fprintf('======================================================\n');
    fprintf('Select an experiment to run:\n');
    fprintf('  [1] 1-D Column Experiment (Surface vs Flux Flip)\n');
    fprintf('  [2] EOF Analysis (Depth x Size Patterns)\n');
    fprintf('  [3] 0-D Slab Test (Sanity Check)\n');
    fprintf('======================================================\n');
    
    choice = input('Enter number [1-3]: ');
    
    switch choice
        case 1
            run_column_experiment();
        case 2
            analyze_model_eofs();
        case 3
            onefile_obs_eps_allin();
        otherwise
            fprintf('Invalid selection. Exiting.\n');
    end
end

% =========================================================================
% === EXPERIMENT 1: 1-D COLUMN (Flux Response) ============================
% =========================================================================
function run_column_experiment
    fprintf('\n=== 1-D COLUMN SIMULATION [AGGRESSIVE PHYSICS] ===\n');

    % 1. Load Data
    if ~isfile('epsilon_daily.mat')
        error('Please ensure epsilon_daily.mat is in the folder.');
    end
    S = load('epsilon_daily.mat');
    [t_days, eps_raw] = grab_data_simple(S);
    
    % Use first 60 days
    mask = t_days <= 60; 
    t_run = t_days(mask);
    eps_run = eps_raw(mask);
    
    fprintf('Running simulation for %.1f days (N=%d steps)...\n', ...
            range(t_run), numel(t_run));

    % 2. Configuration
    cfg = SimulationConfig();
    cfg.use_column = true;
    cfg.z_max      = 200;   
    cfg.dz         = 10;    
    
    % --- AGGRESSIVE PHYSICS (The "Hammer") ---
    cfg.disagg_use_nonlinear = true;
    cfg.disagg_kmax_a        = 0.95; % Break 95% of bins
    cfg.disagg_beta          = 1.0;  % High sensitivity
    cfg.disagg_redistribute_p = 1.5; % Shatter into tiny pieces
    
    cfg.use_NPP     = true;
    cfg.NPP_rate    = 2e-4; 
    cfg.NPP_profile = 'constant';
    
    % 3. Run
    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    
    tic;
    sim.run('tspan', t_run);
    runtime = toc;
    fprintf('Simulation finished in %.2f seconds.\n', runtime);
    
    % 4. Plot using Clean Visualization
    out = sim.exportMinimalOutputs();
    make_clean_plots(out, eps_run);
end

% === CLEAN PLOTTING FUNCTION (Lines Only, No Shading) ===
function make_clean_plots(out, eps_series)
    t = out.t_days;
    flux_200m = out.settling_loss; 
    
    if size(out.N,1) == numel(t), N_mat = out.N.'; else, N_mat = out.N; end

    f = figure('Name','Column Experiment Results','Color','w','Position',[50 50 1200 900]);
    
    % =====================================================================
    % PANEL 1: TURBULENCE (High Visibility)
    % =====================================================================
    ax1 = subplot(3,1,1); 
    
    % Plot data as a black line with dots to ensure spikes are visible
    % 1e-12 ensures the log scale doesn't crash on zeros
    semilogy(t, max(eps_series, 1e-13), 'k.-', 'LineWidth', 1.0, 'MarkerSize', 8);
    hold on;
    
    % Add threshold line
    yline(1e-6, 'r--', 'Threshold (10^{-6})', 'LineWidth', 1.5);
    
    ylabel('\epsilon (W kg^{-1})', 'FontSize', 12, 'FontWeight', 'bold'); 
    title('Forcing: Turbulence', 'FontSize', 14);
    
    % Force limits so spikes don't look like flat lines
    ylim([1e-13, 1e-3]); 
    xlim([t(1) t(end)]);
    grid on; box on;
    
    % =====================================================================
    % PANEL 2: THE FLIP (Flux vs Surface Mass)
    % =====================================================================
    ax2 = subplot(3,1,2);
    yyaxis left
    plot(t, sum(N_mat, 1), 'b-', 'LineWidth', 2.0);
    ylabel('Surface Mass (Blue)', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'YColor', [0 0.45 0.74]); 
    
    yyaxis right
    plot(t, flux_200m, 'r-', 'LineWidth', 2.0);
    ylabel('Export Flux (Red)', 'FontSize', 12, 'FontWeight', 'bold');
    set(gca, 'YColor', [0.85 0.33 0.1]); 
    
    title('The "Flip": Surface Retention vs. Deep Export', 'FontSize', 14);
    xlabel('Time (days)', 'FontSize', 12);
    xlim([t(1) t(end)]); grid on; box on;
    
    % =====================================================================
    % PANEL 3: PSD ANOMALY (High Contrast Heatmap)
    % =====================================================================
    ax3 = subplot(3,1,3);
    if ~isempty(out.D_um)
        N_log = log10(max(N_mat, 1e-20));
        N_mean = mean(N_log, 2); 
        N_anom = N_log - N_mean;
        
        % Use pcolor for smooth look
        h = pcolor(t, log10(out.D_um), N_anom);
        set(h, 'EdgeColor', 'none'); 
        
        set(gca,'YDir','normal'); 
        colormap(jet); 
        caxis([-0.3 0.3]); % High contrast
        
        ylabel('log_{10} Diameter (\mum)', 'FontSize', 12); 
        xlabel('Time (days)', 'FontSize', 12);
        title('Particle Size Anomaly (Red=Gain, Blue=Loss)', 'FontSize', 14);
        
        cb = colorbar; cb.Label.String = 'Log_{10} Anomaly';
        xlim([t(1) t(end)]);
    end
    
    % Link axes so zooming works on all panels
    linkaxes([ax1, ax2, ax3], 'x');
end

% =========================================================================
% === EXPERIMENT 2: EOF ANALYSIS (Model Validation) =======================
% =========================================================================
function analyze_model_eofs
    fprintf('\n=== MODEL EOF ANALYSIS [DEPTH vs SIZE] ===\n');

    % 1. Run Internal Simulation (Same Settings as Exp 1)
    [out, cfg, eps_run] = run_simulation_internal();
    t = out.t_days;
    Y = out.concentrations; 
    
    Nz = floor(cfg.z_max / cfg.dz);
    Ns = cfg.n_sections;
    Nt = length(t);
    
    fprintf('Data: Time=%d, Depth=%d, Size=%d\n', Nt, Nz, Ns);

    % 2. Reshape to [Time, Depth, Size]
    Y_3d = zeros(Nt, Nz, Ns);
    for k = 1:Nz
        idx_start = (k-1)*Ns + 1; idx_end = k*Ns;
        Y_3d(:, k, :) = Y(:, idx_start:idx_end);
    end
    
    % 3. PCA on Log-Anomalies
    Y_flat = reshape(Y_3d, Nt, Nz*Ns); 
    Y_log = log10(max(Y_flat, 1e-20));
    Y_anom = Y_log - mean(Y_log, 1);
    
    [coeffs, scores, ~, ~, explained] = pca(Y_anom);
    
    % 4. Visualization
    d_um = out.D_um; 
    depths = (1:Nz) * cfg.dz; 

    f = figure('Name','EOF Analysis','Color','w','Position',[50 50 1400 900]);
    
    % Panel A: Forcing
    ax1 = subplot(4, 1, 1);
    semilogy(t, max(eps_run,1e-13), 'k.-', 'LineWidth', 0.5);
    yline(1e-6, 'r--'); ylabel('\epsilon'); title('Forcing: Turbulence');
    ylim([1e-13, 1e-3]); axis tight; grid on;
    
    % Panel B: PC Time Series
    ax2 = subplot(4, 1, 2);
    plot(t, normalize(scores(:,2),'range'), 'r-', 'LineWidth', 1.5);
    ylabel('PC2 Amplitude');
    title(sprintf('PC2 (%.1f%%) - The "Flip" Time Series', explained(2)));
    axis tight; grid on;
    
    linkaxes([ax1, ax2], 'x');
    
    % Panel C: EOF 1 (General Mass)
    EOF1_map = reshape(coeffs(:,1), Ns, Nz).'; % [Depth x Size]
    
    subplot(2, 2, 3);
    imagesc(log10(d_um), depths, EOF1_map);
    set(gca, 'YDir', 'reverse'); colormap(jet); colorbar;
    xlabel('log_{10} Diameter'); ylabel('Depth (m)');
    title(sprintf('EOF Mode 1 (%.1f%%)', explained(1)));
    
    % Panel D: EOF 2 (The Signal)
    EOF2_map = reshape(coeffs(:,2), Ns, Nz).'; % [Depth x Size]
    
    subplot(2, 2, 4);
    imagesc(log10(d_um), depths, EOF2_map);
    set(gca, 'YDir', 'reverse'); colormap(jet); colorbar;
    xlabel('log_{10} Diameter'); ylabel('Depth (m)');
    title(sprintf('EOF Mode 2 (%.1f%%) - Disaggregation Signal', explained(2)));
end

function [out, cfg, eps_run] = run_simulation_internal()
    if ~isfile('epsilon_daily.mat'), error('Missing data file'); end
    S = load('epsilon_daily.mat');
    [t, e] = grab_data_simple(S);
    
    mask = t <= 60; % 60 days
    t_run = t(mask);
    eps_run = e(mask);
    
    cfg = SimulationConfig();
    cfg.use_column = true;
    cfg.z_max = 200; cfg.dz = 10;
    
    % Aggressive Physics
    cfg.disagg_use_nonlinear = true;
    cfg.disagg_kmax_a        = 0.95; 
    cfg.disagg_beta          = 1.0; 
    cfg.disagg_redistribute_p = 1.5; 
    cfg.use_NPP = true; cfg.NPP_rate = 2e-4;
    
    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    sim.run('tspan', t_run);
    
    out = sim.exportMinimalOutputs();
    out.concentrations = sim.result.concentrations; 
end

% =========================================================================
% === EXPERIMENT 3: 0-D SLAB TEST (Sanity Check) ==========================
% =========================================================================
function onefile_obs_eps_allin
    fprintf('\n=== 0-D SLAB CHECK [AGGRESSIVE PHYSICS] ===\n');
    if ~isfile('epsilon_daily.mat'), error('File not found'); end
    S = load('epsilon_daily.mat');
    [t_days, eps_raw] = grab_data_simple(S);

    T_total = 30; N1 = 60;
    t1 = linspace(0, T_total, N1).';
    eps1 = interp1(t_days, eps_raw, t1, 'linear', 'extrap');
    
    cfg = SimulationConfig();
    cfg.use_column = false; % FORCE 0-D
    cfg.t_init=0; cfg.t_final=T_total; cfg.delta_t=mean(diff(t1));
    
    % Aggressive Physics
    cfg.disagg_use_nonlinear = true; 
    cfg.disagg_kmax_a = 0.95; cfg.disagg_beta = 1.0; cfg.disagg_redistribute_p = 1.5;
    cfg.use_NPP = true; cfg.NPP_rate = 2e-4;

    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t1, eps1);
    sim.run('tspan', t1);
    out = sim.exportMinimalOutputs();
    
    if isempty(out.D_um)
        out.D_um = 2 * sim.grid.getFractalRadii() * 1e4; 
    end
    
    f = figure('Name','Slab Check','Color','w');
    subplot(2,1,1); semilogy(t1, max(eps1, 1e-12), 'k.-'); title('Turbulence'); ylim([1e-13 1e-3]);
    
    N_log = log10(max(out.N.', 1e-20));
    N_anom = N_log - mean(N_log, 2);
    
    subplot(2,1,2); imagesc(t1, log10(out.D_um), N_anom);
    set(gca,'YDir','normal'); colormap(jet); caxis([-0.5 0.5]);
    title('Slab Anomaly (Red=Gain, Blue=Loss)'); colorbar;
end

% =========================================================================
% === SHARED HELPERS ======================================================
% =========================================================================
function [t, e] = grab_data_simple(S)
    names = fieldnames(S);
    if numel(names)==1 && isstruct(S.(names{1})), S=S.(names{1}); names=fieldnames(S); end
    t=[]; e=[];
    for k=1:numel(names)
        val = S.(names{k}); nm = lower(names{k});
        if (contains(nm,'time')||contains(nm,'day')) && isnumeric(val), t=val; end
        if (contains(nm,'eps')||contains(nm,'dissip')) && isnumeric(val), e=val; end
    end
    t = t(:) - t(1);
    n = min(numel(t), numel(e));
    t = t(1:n); e = e(1:n);
    e = max(e, 1e-12);
end