function analyze_model_eofs
    fprintf('=== MODEL EOF ANALYSIS [DEPTH vs SIZE] ===\n');

    % 1. Run Simulation
    [out, cfg, eps_run] = run_simulation_internal();
    t = out.t_days;
    Y = out.concentrations; 
    
    Nz = floor(cfg.z_max / cfg.dz);
    Ns = cfg.n_sections;
    Nt = length(t);
    
    fprintf('Data: Time=%d, Depth=%d, Size=%d\n', Nt, Nz, Ns);

    % 2. Reshape [Time, Depth, Size]
    Y_3d = zeros(Nt, Nz, Ns);
    for k = 1:Nz
        idx_start = (k-1)*Ns + 1;
        idx_end   = k*Ns;
        Y_3d(:, k, :) = Y(:, idx_start:idx_end);
    end
    
    % 3. PCA on Log-Anomalies
    Y_flat = reshape(Y_3d, Nt, Nz*Ns); 
    Y_log = log10(max(Y_flat, 1e-20));
    Y_mean = mean(Y_log, 1);
    Y_anom = Y_log - Y_mean;
    
    [coeffs, scores, ~, ~, explained] = pca(Y_anom);
    
    % 4. Recover Diameters for Axis
    d_um = out.D_um; 
    depths = (1:Nz) * cfg.dz; % Depth vector in meters

    % 5. Visualization (Standard Oceanography Layout)
    f = figure('Color','w','Position',[50 50 1400 900]);
    
    % --- Panel A: Forcing ---
    subplot(4, 1, 1);
    semilogy(t, eps_run, 'k-', 'LineWidth', 1.2);
    yline(1e-6, 'r--', 'Threshold');
    ylabel('\epsilon (W/kg)'); title('Forcing: Turbulence');
    axis tight; grid on;
    
    % --- Panel B: PC Time Series ---
    subplot(4, 1, 2);
    plot(t, normalize(scores(:,2),'range'), 'r-', 'LineWidth', 1.5);
    ylabel('PC2 Amplitude');
    title(sprintf('PC2 (%.1f%%) - The "Flip" Time Series', explained(2)));
    axis tight; grid on;
    
    % --- Panel C & D: EOF Spatial Patterns (Depth vs Size) ---
    % We want: Y-axis = Depth (reversed), X-axis = Size
    
    % EOF 1
    EOF1_flat = coeffs(:,1);
    EOF1_map  = reshape(EOF1_flat, Nz, Ns); % [Depth x Size]
    
    subplot(2, 2, 3);
    imagesc(log10(d_um), depths, EOF1_map);
    set(gca, 'YDir', 'reverse'); % Depth goes down
    colormap(jet); colorbar;
    xlabel('log10 Diameter (\mum)'); ylabel('Depth (m)');
    title(sprintf('EOF Mode 1 (%.1f%%)', explained(1)));
    
    % EOF 2 (The Flip)
    EOF2_flat = coeffs(:,2);
    EOF2_map  = reshape(EOF2_flat, Nz, Ns);
    
    subplot(2, 2, 4);
    imagesc(log10(d_um), depths, EOF2_map);
    set(gca, 'YDir', 'reverse'); % Depth goes down
    colormap(jet); colorbar;
    xlabel('log10 Diameter (\mum)'); ylabel('Depth (m)');
    title(sprintf('EOF Mode 2 (%.1f%%) - The Disaggregation Signal', explained(2)));
    
    fprintf('Done. Check EOF Mode 2 for the vertical dipole pattern.\n');
end

function [out, cfg, eps_run] = run_simulation_internal()
    if ~isfile('epsilon_daily.mat'), error('Missing data file'); end
    S = load('epsilon_daily.mat');
    names = fieldnames(S);
    if numel(names)==1 && isstruct(S.(names{1})), S=S.(names{1}); names=fieldnames(S); end
    t=[]; e=[];
    for k=1:numel(names)
        v=S.(names{k}); nm=lower(names{k});
        if contains(nm,'eps')||contains(nm,'dis'), e=v; end
        if contains(nm,'time')||contains(nm,'day'), t=v; end
    end
    t = t(:) - t(1);
    mask = t <= 60;
    t_run = t(mask);
    eps_run = e(mask);
    
    cfg = SimulationConfig();
    cfg.use_column = true;
    cfg.z_max = 200; cfg.dz = 10;
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