function master_simulation_runner
    % MASTER_SIMULATION_RUNNER
    % All-in-one script for Particle Aggregation Project.
    % Version: square-wave test, flux fractions, EOF + 3D slice plots.
    
    clc;
    fprintf('======================================================\n');
    fprintf('   PARTICLE AGGREGATION MODEL - RESEARCH EDITION\n');
    fprintf('======================================================\n');
    fprintf('Select an experiment:\n');
    fprintf('  [1] Real Atlantic Data (Full Analysis: Relative Flux + Size Classes)\n');
    fprintf('  [2] Square Wave Test (Idealized Lab Test)\n');
    fprintf('  [3] EOF Analysis (Atlantic Validation + 3D Slice Plots)\n');
    fprintf('  [4] 0-D Slab Test (Sanity Check)\n');
    fprintf('  [5] North Pacific Flux Analysis (Compare to Atlantic Flux)\n');
    fprintf('======================================================\n');
    
    choice = input('Enter number [1-5]: ');
    
    switch choice
        case 1
            run_column_experiment();
        case 2
            run_square_wave_experiment();
        case 3
            analyze_model_eofs('Atlantic');
        case 4
            onefile_obs_eps_allin();
        case 5 
            run_north_pacific_experiment(); 
        otherwise
            fprintf('Invalid selection.\n');
    end
end

% =========================================================================
% === EXPERIMENT 1: REAL ATLANTIC DATA ====================================
% =========================================================================
function run_column_experiment
    fprintf('\n=== 1-D COLUMN SIMULATION [REAL FORCING + ATTENUATION] ===\n');
    if ~isfile('epsilon_daily.mat')
        error('Please ensure epsilon_daily.mat is in the folder.');
    end
    S = load('epsilon_daily.mat');
    [t_days, eps_raw] = grab_data_simple(S);
    
    mask    = t_days <= 60; 
    t_run   = t_days(mask);
    eps_run = eps_raw(mask);
    
    fprintf('Running simulation for %.1f days (N=%d steps)...\n', ...
        range(t_run), numel(t_run));
    cfg = SimulationConfig();
    cfg.use_column = true; cfg.z_max = 200; cfg.dz = 10;
    cfg.disagg_use_nonlinear = true; 
    cfg.disagg_kmax_a = 0.95; 
    cfg.disagg_beta   = 1.0; 
    cfg.disagg_redistribute_p = 1.5;
    cfg.attenuation_rate = 0.1; 
    cfg.use_NPP = true; cfg.NPP_rate = 2e-4;
    
    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    
    tic; sim.run('tspan', t_run); 
    fprintf('Simulation finished in %.2f seconds.\n', toc);
    
    make_comprehensive_plots(sim.result, sim.grid, cfg, eps_run, ...
        'Atlantic Observed');
end

% =========================================================================
% === EXPERIMENT 2: SQUARE WAVE TEST ======================================
% =========================================================================
function run_square_wave_experiment
    fprintf('\n=== SQUARE WAVE EXPERIMENT (WITH ATTENUATION) ===\n');
    
    T_total = 15;
    dt      = 0.05;              % fine time resolution
    t_run   = (0:dt:T_total).';
    
    % Background low epsilon
    eps_run = 1e-9 * ones(size(t_run));
    
    % Three square pulses (~1 day wide) like Adrian's sketch
    pulse_edges = [0.5 1.5;
                   5.0 6.0;
                   9.5 10.5];
    for k = 1:size(pulse_edges,1)
        mask = (t_run >= pulse_edges(k,1)) & (t_run < pulse_edges(k,2));
        eps_run(mask) = 5e-6;
    end
    
    cfg = SimulationConfig();
    cfg.use_column = true; cfg.z_max = 200; cfg.dz = 10;
    cfg.disagg_use_nonlinear = true; 
    cfg.disagg_kmax_a = 0.95; 
    cfg.disagg_beta   = 1.0; 
    cfg.disagg_redistribute_p = 1.5;
    cfg.attenuation_rate = 0.02; % slightly weaker for clearer overshoot
    cfg.use_NPP = true; cfg.NPP_rate = 2e-4;
    
    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    sim.run('tspan', t_run);
    
    make_comprehensive_plots(sim.result, sim.grid, cfg, eps_run, ...
        'Square Wave Test');
end

% =========================================================================
% === EXPERIMENT 5: NORTH PACIFIC FLUX ANALYSIS ===========================
% =========================================================================
function run_north_pacific_experiment
    fprintf('\n=== NORTH PACIFIC FLUX ANALYSIS (BIOLOGY-DRIVEN) ===\n');
    data_file = 'epsilon_NP.mat'; 
    if ~isfile(data_file)
        error('Please ensure the North Pacific data file "%s" is in the folder.', ...
            data_file);
    end
    S = load(data_file);
    [t_days, eps_raw] = grab_data_simple(S);
    
    mask    = t_days <= 60; 
    t_run   = t_days(mask);
    eps_run = eps_raw(mask);
    
    fprintf('Running simulation for %.1f days (N=%d steps)...\n', ...
        range(t_run), numel(t_run));
    cfg = SimulationConfig();
    cfg.use_column = true; cfg.z_max = 200; cfg.dz = 10;
    cfg.disagg_use_nonlinear = true; 
    cfg.disagg_kmax_a = 0.95; 
    cfg.disagg_beta   = 1.0; 
    cfg.disagg_redistribute_p = 1.5;
    cfg.attenuation_rate = 0.1; 
    cfg.use_NPP = true; cfg.NPP_rate = 2e-4;
    
    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    
    tic; sim.run('tspan', t_run); 
    fprintf('Simulation finished in %.2f seconds.\n', toc);
    
    make_comprehensive_plots(sim.result, sim.grid, cfg, eps_run, ...
        'North Pacific');
end

% =========================================================================
% === CORE PLOTTING FUNCTION (FLUX & SIZE CLASSES) ========================
% =========================================================================
function make_comprehensive_plots(result, sim_grid_obj, cfg, eps_run, ...
    title_prefix)

    t  = result.time; 
    Y  = result.concentrations;
    Nt = length(t); 
    Nz = floor(cfg.z_max / cfg.dz); 
    Ns = cfg.n_sections;
    
    r_cm = sim_grid_obj.getFractalRadii(); 
    r_v  = sim_grid_obj.getConservedRadii();
    ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, sim_grid_obj.setcon);
    W_m_d   = (ws_cm_s / 100) * 86400;         % cm/s -> m/day
    V_m3    = (4/3) * pi * r_v.^3 * 1e-6;      % cm^3 -> m^3
    D_um    = 2 * r_cm * 1e4;
    
    export_depth_idx = Nz; 
    Q_export = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);
    F_total  = sum(Q_export .* V_m3.' .* W_m_d.', 2); 
    
    % Steady-state flux (background epsilon)
    F_ss_idx = find(eps_run < 1e-10, 10); 
    if isempty(F_ss_idx), F_ss_idx = 1:min(10,Nt); end 
    F_ss = mean(F_total(F_ss_idx)); 
    if F_ss == 0, F_ss = 1; end 
    F_relative = F_total / F_ss;
    
    % Size classes (micron)
    D_small_max  = 200; 
    D_medium_max = 1000; 
    idx_small  = (D_um < D_small_max);
    idx_medium = (D_um >= D_small_max) & (D_um < D_medium_max);
    idx_large  = (D_um >= D_medium_max);
    
    F_small  = sum(Q_export(:, idx_small)  .* V_m3(idx_small).'  .* ...
                   W_m_d(idx_small).',  2);
    F_medium = sum(Q_export(:, idx_medium) .* V_m3(idx_medium).' .* ...
                   W_m_d(idx_medium).', 2);
    F_large  = sum(Q_export(:, idx_large)  .* V_m3(idx_large).'  .* ...
                   W_m_d(idx_large).',  2);
    
    % Safe total for fractions
    F_total_safe = max(F_total, 1e-30);
    frac_small   = F_small  ./ F_total_safe;
    frac_medium  = F_medium ./ F_total_safe;
    frac_large   = F_large  ./ F_total_safe;
    
    figure('Name', [title_prefix ' Flux Analysis'], ...
           'Color', 'w', 'Position', [100 100 800 900]);
    
    % Panel 1: turbulence
    subplot(3, 1, 1);
    semilogy(t, max(eps_run, 1e-12), 'k-'); hold on;
    plot(t([1 end]), [1e-6 1e-6], 'r--', 'LineWidth', 1); 
    ylabel('\epsilon (W kg^{-1})');
    title([title_prefix ' Simulation']); 
    axis tight; grid on;
    
    % Panel 2: relative total export flux
    subplot(3, 1, 2);
    plot(t, F_relative, 'b-', 'LineWidth', 2); hold on;
    plot(t([1 end]), [1 1], 'k--'); 
    ylabel('Relative Flux (F/F_{ss})');
    title('Relative Export Flux');
    axis tight; grid on;
    
    % Panel 3: flux fractions (Large/Medium/Small) – linear axis
    subplot(3, 1, 3);
    plot(t, frac_large,  'r-',  'LineWidth', 2.2, ...
             'DisplayName', 'Large');
    hold on;
    plot(t, frac_medium, 'g--', 'LineWidth', 2.2, ...
             'DisplayName', 'Medium');
    plot(t, frac_small,  'b:',  'LineWidth', 2.2, ...
             'DisplayName', 'Small');
    ylim([0 1]);
    ylabel('Flux Fraction');
    xlabel('Time (days)');
    legend('Location', 'southwest');
    grid on; box on;
end

% =========================================================================
% === UNIFIED EOF ANALYSIS + 3D SLICE PLOTS ===============================
% =========================================================================
function analyze_model_eofs(title_prefix, t, Y, cfg, eps_run)
    % If called with only the title, run an internal Atlantic simulation
    if nargin < 2
         fprintf('\n=== MODEL EOF ANALYSIS [WITH ATTENUATION] ===\n');
         [out, cfg, eps_run] = run_simulation_internal();
         t            = out.t_days; 
         Y            = out.concentrations; 
         title_prefix = 'Atlantic';
    end
    
    Nz = floor(cfg.z_max / cfg.dz); 
    Ns = cfg.n_sections; 
    Nt = length(t);
    depths = (1:Nz) * cfg.dz; 

    % --- reshape to [time, depth, size] ---------------------------------
    Y_3d = zeros(Nt, Nz, Ns);
    for k = 1:Nz
        idx_start = (k-1)*Ns + 1; 
        idx_end   = k*Ns;
        Y_3d(:, k, :) = Y(:, idx_start:idx_end);
    end
    
    % --- log10 PSD, anomalies for EOFs ----------------------------------
    Y_flat = reshape(Y_3d, Nt, Nz*Ns); 
    Y_log  = log10(max(Y_flat, 1e-20));   % log10(PartVol or number)
    Y_mean = mean(Y_log, 1);
    Y_anom = Y_log - Y_mean;
    
    [coeffs, scores, ~, ~, explained] = pca(Y_anom);   % coeffs: [Nz*Ns x nModes]
    
    % --- size axis in mm -------------------------------------------------
    if exist('out','var')
        D_um = out.D_um;
    else
        sim_temp = CoagulationSimulation(cfg); 
        D_um     = 2 * sim_temp.grid.getFractalRadii() * 1e4;
    end
    D_mm  = D_um / 1000;
    logD  = log10(D_mm);
    
    % --- reshape EOF modes into [depth x size] ---------------------------
    mode1 = reshape(coeffs(:,1), Ns, Nz).';   % [Nz x Ns]
    mode2 = reshape(coeffs(:,2), Ns, Nz).';
    
    % Flip sign of Mode 1 so the surface / large sizes are positive
    if mean(mode1(1, :)) < 0
        mode1       = -mode1;
        coeffs(:,1) = -coeffs(:,1);
        scores(:,1) = -scores(:,1);
    end
    
    % Flip Mode 2 so that large particles near the surface are positive
    if mean(mode2(1, end-2:end)) < 0
        mode2       = -mode2;
        coeffs(:,2) = -coeffs(:,2);
        scores(:,2) = -scores(:,2);
    end
    
    % --- FIGURE 1: EOF maps + forcing / PCs ------------------------------
    figure('Name',[title_prefix ' EOF Analysis'],'Color','w', ...
           'Position',[50 50 1200 850]);
    
    % (a) Turbulence forcing
    subplot(2, 2, 1); 
    semilogy(t, max(eps_run,1e-13), 'k-','LineWidth',1.2); 
    ylabel('\epsilon (W kg^{-1})'); 
    title('Turbulence Forcing'); 
    axis tight; grid on;
    
    % (b) PC1 & PC2 time series (normalized)
    subplot(2, 2, 2); 
    pc1 = scores(:,1);
    pc2 = scores(:,2);
    pc1_n = pc1 ./ max(abs(pc1));
    pc2_n = pc2 ./ max(abs(pc2));
    plot(t, pc1_n, 'k-', 'LineWidth', 1.2, 'DisplayName','PC1');
    hold on;
    plot(t, pc2_n, 'r-', 'LineWidth', 1.2, 'DisplayName','PC2');
    yline(0,'k--');
    ylabel('PC (normalized)');
    title(sprintf('PC1 & PC2 (%.1f%% / %.1f%% var)', explained(1), explained(2)));
    legend('Location','best');
    axis tight; grid on;
    
    % (c) Mode 1 map (attenuation)
    subplot(2, 2, 3);
    imagesc(logD, depths, mode1); 
    set(gca,'YDir','reverse');
    colormap(gca, turbo); 
    cb = colorbar;
    cb.Label.String = 'EOF1 loading';
    title(sprintf('Mode 1 (Attenuation, %.1f%%%% var)', explained(1)));
    xlabel('D (mm)'); ylabel('Depth (m)');
    set(gca,'XTick',log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});
    
    % (d) Mode 2 map (seesaw)
    subplot(2, 2, 4);
    imagesc(logD, depths, mode2); 
    set(gca,'YDir','reverse');
    colormap(gca, turbo); 
    cb = colorbar; 
    cb.Label.String = 'EOF2 loading';
    title(sprintf('Mode 2 (Seesaw, %.1f%%%% var)', explained(2)));
    xlabel('D (mm)'); ylabel('Depth (m)');
    set(gca,'XTick',log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});
    
    % --- FIGURE 2: 3D slice plots (Full, Mode1, Mode2) -------------------
    PSD_full_log = Y_log;  
    PSD_full_3d  = reshape(PSD_full_log, Nt, Nz, Ns);
    
    Y1_anom = scores(:,1) * coeffs(:,1).';      
    Y1_log  = Y1_anom + Y_mean;                 
    Y1_3d   = reshape(Y1_log, Nt, Nz, Ns);
    
    Y2_anom = scores(:,2) * coeffs(:,2).';
    Y2_log  = Y2_anom + Y_mean; 
    Y2_3d   = reshape(Y2_log, Nt, Nz, Ns);

    % colour limits from 2nd–98th percentiles (so we see structure)
    vals = sort(PSD_full_log(:));
    vals = vals(isfinite(vals));
    n    = numel(vals);
    c_min = vals(max(1, round(0.02*n)));
    c_max = vals(round(0.98*n));
    clim_global = [c_min c_max];

    figure('Name', [title_prefix ' 3D Reconstruction'], 'Color', 'w', ...
           'Position', [50 50 1100 950]);
    
    ax1 = subplot(3,1,1);
    plot_psd_slice_subplot(ax1, t, depths, D_mm, PSD_full_3d, ...
        'Full  Model (log_{10} Part Vol)', clim_global);
    
    ax2 = subplot(3,1,2);
    plot_psd_slice_subplot(ax2, t, depths, D_mm, Y1_3d, ...
        'Mode 1', clim_global);
    
    ax3 = subplot(3,1,3);
    plot_psd_slice_subplot(ax3, t, depths, D_mm, Y2_3d, ...
        'Mode 2', clim_global);
end

% =========================================================================
% === HELPER: 3D SLICE PLOT (FENCES) – SIDE COLORBAR, TIME ORDER ==========
% =========================================================================
function plot_psd_slice_subplot(ax_handle, t_days, depths, D_mm, PSD_3d, ...
    title_str, c_limits)

    axes(ax_handle); 
    hold on;

    t_days = t_days(:);
    depths = depths(:);
    logD   = log10(D_mm(:));      % internal coord (log10 mm)

    % PSD_3d is [Nt, Nz, Ns] -> V [Ns, Nt, Nz] for slice()
    V = permute(PSD_3d, [3, 1, 2]); 
    [X, Y, Z] = meshgrid(t_days, logD, depths); 

    % one "fence" per day; flip so earliest days are at the front
    t_min  = ceil(min(t_days));
    t_max  = floor(max(t_days));
    xs = t_min : 1 : t_max;       % change to :2: for fewer fences
    xslice = fliplr(xs);
    yslice = [];
    zslice = [];

    h = slice(X, Y, Z, V, xslice, yslice, zslice);
    set(h, 'EdgeColor', 'none', 'FaceAlpha', 0.95);

    set(gca, 'ZDir', 'reverse'); 
    colormap(gca, turbo); 
    caxis(c_limits); 
    
    % vertical colourbar on the side (compact)
    cb = colorbar('Location','eastoutside');
    cb.Label.String = 'log_{10} Part Vol (ppmV mm^{-1})';

    xlabel('Time (days)');
    ylabel('D (mm)');
    zlabel('Depth (m)');
    title(title_str);
    
    % log D axis with ticks at 0.1, 1, 10 mm
    set(gca,'YTick', log10([0.1 1 10]), ...
            'YTickLabel', {'0.1','1','10'});
    
    % stretch axis so fences are clear, orientation like data figure
    pbaspect([35 1 3]);          % [x(time) y(D) z(depth)]
    view([-60 5]);               % earliest days in front
    axis tight; box on; grid on;
end

% === UTILITY FUNCTIONS ===================================================
function [out, cfg, eps_run] = run_simulation_internal()
    if ~isfile('epsilon_daily.mat')
        error('Missing epsilon_daily.mat for internal EOF run.');
    end
    S = load('epsilon_daily.mat'); 
    [t, e] = grab_data_simple(S);
    mask    = t <= 60; 
    t_run   = t(mask); 
    eps_run = e(mask);
    
    cfg = SimulationConfig();
    cfg.use_column = true; cfg.z_max = 200; cfg.dz = 10;
    cfg.disagg_use_nonlinear = true; 
    cfg.disagg_kmax_a = 0.95; 
    cfg.disagg_beta   = 1.0; 
    cfg.disagg_redistribute_p = 1.5;
    cfg.attenuation_rate = 0.01; 
    cfg.use_NPP = true; cfg.NPP_rate = 2e-4;
    
    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    sim.run('tspan', t_run);
    
    out_temp           = sim.exportMinimalOutputs();
    out.t_days         = sim.result.time;
    out.concentrations = sim.result.concentrations; 
    out.D_um           = out_temp.D_um;
end

function onefile_obs_eps_allin
    fprintf('Skipping 0-D check for now.\n');
end

function [t, e] = grab_data_simple(S)
    names = fieldnames(S);
    if numel(names)==1 && isstruct(S.(names{1}))
        S     = S.(names{1}); 
        names = fieldnames(S); 
    end
    t = []; e = [];
    for k = 1:numel(names)
        val = S.(names{k}); 
        nm  = lower(names{k});
        if (contains(nm,'time') || contains(nm,'day')) && isnumeric(val)
            t = val; 
        end
        if (contains(nm,'eps') || contains(nm,'dissip')) && isnumeric(val)
            e = val; 
        end
    end
    t = t(:) - t(1);
    n = min(numel(t), numel(e));
    t = t(1:n); 
    e = e(1:n);
    e = max(e, 1e-12);   % floor to avoid zeros
end