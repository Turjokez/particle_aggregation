function run_column_experiment
    fprintf('=== 1-D COLUMN SIMULATION: OBSERVED TURBULENCE [AGGRESSIVE] ===\n');

    % 1. Load Data
    if ~isfile('epsilon_daily.mat')
        error('Please ensure epsilon_daily.mat is in the folder.');
    end
    S = load('epsilon_daily.mat');
    [t_days, eps_raw] = grab_data_simple(S);
    
    % Use first 60 days to see the first major storms clearly
    mask = t_days <= 60; 
    t_run = t_days(mask);
    eps_run = eps_raw(mask);
    
    fprintf('Running simulation for %.1f days (N=%d steps)...\n', ...
            range(t_run), numel(t_run));

    % 2. Configuration
    cfg = SimulationConfig();
    
    % --- 1-D COLUMN SETTINGS ---
    cfg.use_column = true;
    cfg.z_max      = 200;   % 200 meters deep
    cfg.dz         = 10;    % 10m layers
    
    % --- PHYSICS: AGGRESSIVE BREAKUP (The "Hammer") ---
    cfg.disagg_use_nonlinear = true;
    
    % Break 95% of large size classes when turbulence is high
    cfg.disagg_kmax_a        = 0.95; 
    
    % Reaction speed: Very sensitive to epsilon spikes
    cfg.disagg_beta          = 1.0; 
    
    % Shatter into DUST: Force mass into very small bins
    cfg.disagg_redistribute_p = 1.5; 
    
    % --- BIOLOGY ---
    cfg.use_NPP     = true;
    cfg.NPP_rate    = 2e-4; 
    cfg.NPP_profile = 'constant';
    
    % 3. Initialize & Attach Data
    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    
    % 4. Run
    tic;
    sim.run('tspan', t_run);
    runtime = toc;
    fprintf('Simulation finished in %.2f seconds.\n', runtime);
    
    % 5. Extract Results
    out = sim.exportMinimalOutputs();
    
    % 6. Plotting
    make_column_plots(out, eps_run);
end

% =========================================================================
% PLOTTING FUNCTION
% =========================================================================
function make_column_plots(out, eps_series)
    t = out.t_days;
    flux_200m = out.settling_loss; 
    
    if size(out.N,1) == numel(t)
        mass_surface = sum(out.N, 2);
        N_mat = out.N.'; 
    else
        mass_surface = sum(out.N, 1).';
        N_mat = out.N;
    end

    f = figure('Color','w','Position',[50 50 1000 700]);
    
    % Panel 1: Turbulence Input
    subplot(3,1,1);
    semilogy(t, eps_series, 'k-');
    yline(1e-6,'r--','Threshold');
    ylabel('\epsilon (W/kg)');
    title('Forcing: Turbulence');
    grid on; axis tight;
    
    % Panel 2: The "Flip" Check
    subplot(3,1,2);
    yyaxis left
    plot(t, mass_surface, 'b-', 'LineWidth', 1.5);
    ylabel('Surface Mass (Blue)', 'Color', 'b');
    
    yyaxis right
    plot(t, flux_200m, 'r-', 'LineWidth', 1.5);
    ylabel('Export Flux (Red)', 'Color', 'r');
    
    title('The "Flip": Does Flux Drop when Surface Mass Spikes?');
    xlabel('Time (days)');
    grid on; axis tight;
    
    % Panel 3: Surface PSD Heatmap
    subplot(3,1,3);
    if ~isempty(out.D_um)
        imagesc(t, log10(out.D_um), log10(max(N_mat, 1e-20)));
        set(gca,'YDir','normal');
        colormap(jet); 
        
        % === VISUALIZATION FIX ===
        % Tighten colors to show the "gap" in large particles
        caxis([-3.5 -0.5]); 
        
        ylabel('log10 Diameter (\mum)');
        xlabel('Time (days)');
        title('Surface PSD (Look for vertical stripes/gaps)');
        colorbar;
    end
end

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
end