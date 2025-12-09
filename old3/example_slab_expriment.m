%% === example_slab_experiment.m ===

% --------------------------------------------------------------
clear; clc; close all;
fprintf('\n=== Example Slab Experiment — Full Combined Version ===\n');

outdir = fullfile(pwd, 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end
save_figs = true;
colors = lines(6);
lw = 1.2; fs = 11;

%% === BASE CONFIGURATION ===
base_cfg = SimulationConfig();
base_cfg.t_final     = 30;
base_cfg.use_NPP     = true;
base_cfg.NPP_rate    = 5e-4;
base_cfg.NPP_profile = 'sine'; % Sinusoidal NPP forcing (±50%, 30-day period)
base_cfg.epsilon     = 1e-6;
base_cfg.temp = 20; 
base_cfg.salt = 35;

fprintf('Configuration ready → running baseline simulation...\n');
cfg  = base_cfg;                          % copy for baseline
sim  = CoagulationSimulation(cfg);
out  = sim.run();
grid = DerivedGrid(cfg);

fprintf('Simulation finished: %d time steps × %d sections\n',...
    numel(out.time), size(out.concentrations,2));

%% === FIG 1: Baseline (NPP, ε, Flux) ===
tt = linspace(0, cfg.t_final, 100);
NPP_t = cfg.NPP_rate * (1 + 0.5*sin(2*pi*tt/cfg.t_final));
eps_t = cfg.epsilon * ones(size(tt));
[t_flux, flux] = OutputGenerator.simpleFluxTimeSeries(out, grid);

figure('Name','Core_Baseline','Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
nexttile; plot(tt,NPP_t,'b','LineWidth',lw); grid on;
xlabel('Time (d)'); ylabel('NPP (model d^{-1})');
title('NPP vs Time','FontWeight','bold');
nexttile; semilogy(tt,eps_t,'r','LineWidth',lw); grid on;
xlabel('Time (d)'); ylabel('\epsilon (W kg^{-1})');
title('\epsilon vs Time','FontWeight','bold'); ylim([1e-8 1e-5]);
nexttile; plot(out.time,flux,'k','LineWidth',lw); grid on;
xlabel('Time (d)'); ylabel('Flux (a.u.)');
title('Flux vs Time','FontWeight','bold');
if save_figs, exportgraphics(gcf,fullfile(outdir,'Core_Baseline.png'),'Resolution',300); end


%% === FIG 2: ε SWEEP — Flux vs Time + τ50 ===
fprintf('\n[1] Running ε sweep...\n');
eps_list = [1e-8 1e-7 1e-6 1e-5];
tau50_rise = @(t,f) compute_tau50_rise(t,f);
tau50_eps = NaN(size(eps_list));
flux_curves = cell(numel(eps_list),1);
t_curves = cell(numel(eps_list),1);
Y_final = cell(numel(eps_list),1);
grid_all = cell(numel(eps_list),1);

for k = 1:numel(eps_list)
    cfgk = base_cfg; 
    cfgk.epsilon = eps_list(k);             % only vary ε
    simk = CoagulationSimulation(cfgk);
    outk = simk.run();
    gridk = DerivedGrid(cfgk);
    [tk,fk] = OutputGenerator.simpleFluxTimeSeries(outk,gridk);
    flux_curves{k} = fk; 
    t_curves{k} = tk;
    tau50_eps(k) = tau50_rise(tk,fk);
    Y_final{k} = outk.concentrations(end,:)'; 
    grid_all{k} = gridk;
end

figure('Name','eps_sweep','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');
nexttile; hold on;
for k=1:numel(eps_list)
    plot(t_curves{k}, flux_curves{k}, 'LineWidth',lw, 'Color',colors(k,:));
end
hold off; grid on;
xlabel('Time (d)'); ylabel('Flux (a.u.)');
title('Flux vs Time — \epsilon sweep','FontWeight','bold');
legend(arrayfun(@(e)sprintf('\\epsilon=%.0e',e),eps_list,'UniformOutput',false),'Location','best');
nexttile; mask=isfinite(tau50_eps);
if any(mask)
    semilogx(eps_list(mask), tau50_eps(mask),'o-k','LineWidth',lw,'MarkerFaceColor','w');
else
    text(1e-7,0.5,'No valid \tau_{50}','FontSize',fs);
end
xlabel('\epsilon (W kg^{-1})'); ylabel('\tau_{50} (d)');
title('\tau_{50} vs \epsilon','FontWeight','bold'); grid on;
if save_figs, exportgraphics(gcf,fullfile(outdir,'eps_sweep.png'),'Resolution',300); end


%% === FIG 3: Final PSDs for ε cases ===
figure('Name','epsilon_sweep_PSD','Color','w'); hold on;
nE = numel(eps_list);
slopes = nan(1,nE);

for k = 1:nE
    if k > numel(Y_final) || isempty(Y_final{k})
        continue;
    end
    gridk = grid_all{k};
    diam = 2 * gridk.getFractalRadii() * 1e4;   % cm → µm
    y = Y_final{k};
    valid = y > 0 & isfinite(y);
    if nnz(valid) < 3, slopes(k)=NaN; continue; end
    D = diam(valid); N = y(valid);
    p = polyfit(log10(D), log10(N), 1);
    b = -p(1);
    slopes(k) = b;
    fit_line = 10.^(polyval(p, log10(D)));
    loglog(D, N, 'LineWidth', 1.5, 'Color', colors(k,:));
    loglog(D, fit_line, '--', 'LineWidth', 1.1, 'Color', colors(k,:)*0.7);
end

xlabel('Particle diameter (µm)','FontSize',11);
ylabel('Concentration (a.u.)','FontSize',11);
title('Final Particle Size Distribution vs Turbulence (\epsilon)','FontWeight','bold');

mask = isfinite(slopes);
legend_entries = arrayfun(@(x,b) sprintf('\\epsilon = %.0e  (B = %.2f)', x, b), ...
    eps_list(mask), slopes(mask), 'UniformOutput', false);
legend(legend_entries, 'Location','southwest','FontSize',10);
grid on; box on;
set(gca,'XScale','log','YScale','log',...
    'XLim',[1e0 1e4],'YLim',[1e-10 1e-4],...
    'LineWidth',1,'FontSize',10,'MinorGridLineStyle',':');

exportgraphics(gcf, fullfile(outdir,'epsilon_sweep_PSD_fitted.png'),'Resolution',300);
fprintf('\n✅ PSD figure with fitted slopes saved to: %s\n', fullfile(outdir,'epsilon_sweep_PSD_fitted.png'));


%% === FIG 4: NPP RESPONSE TEST ===
fprintf('\n[2] Running NPP sensitivity test...\n');
npp_levels = [1e-4 5e-4 1e-3];
tN = cell(numel(npp_levels),1); 
fN = cell(numel(npp_levels),1);
for k=1:numel(npp_levels)
    cfgk = base_cfg; 
    cfgk.NPP_rate    = npp_levels(k); 
    cfgk.NPP_profile = 'constant'; 
    cfgk.epsilon     = base_cfg.epsilon;  % fixed turbulence
    outk = CoagulationSimulation(cfgk).run();
    [tN{k}, fN{k}] = OutputGenerator.simpleFluxTimeSeries(outk, DerivedGrid(cfgk));
end

figure('Name','npp_response','Color','w'); hold on;
for k=1:numel(npp_levels)
    plot(tN{k}, fN{k}, 'LineWidth', lw, 'Color', colors(k,:));
end
xlabel('Time (d)'); ylabel('Flux (a.u.)');
title('Flux vs Time — NPP levels','FontWeight','bold');
legend(arrayfun(@(x)sprintf('NPP=%.1e',x),npp_levels,'UniformOutput',false),'Location','best');
grid on;
if save_figs, exportgraphics(gcf,fullfile(outdir,'npp_response.png'),'Resolution',300); end


%% === FIG 5: MLD PROXY TEST ===
fprintf('\n[3] Running MLD proxy test (10–50 m)...\n');
mld_list = [10 20 30 40 50]; 
ref_mld  = 30;
tM = cell(numel(mld_list),1); 
fM = cell(numel(mld_list),1);

for k=1:numel(mld_list)
    cfgk = base_cfg;
    scale = ref_mld / mld_list(k);
    cfgk.NPP_profile = 'constant';
    cfgk.NPP_rate    = base_cfg.NPP_rate * scale;
    cfgk.epsilon     = base_cfg.epsilon;
    outk = CoagulationSimulation(cfgk).run();
    [tM{k}, fM{k}] = OutputGenerator.simpleFluxTimeSeries(outk, DerivedGrid(cfgk));
end

figure('Name','mld_proxy','Color','w'); hold on;
for k=1:numel(mld_list)
    plot(tM{k}, fM{k}, 'LineWidth', lw, 'Color', colors(k,:));
end
xlabel('Time (d)'); ylabel('Flux (a.u.)');
title('Flux vs Time — mixed-layer proxy (10–50 m)','FontWeight','bold');
legend(arrayfun(@(z)sprintf('MLD %2.0f m',z),mld_list,'UniformOutput',false),'Location','best');
grid on;
if save_figs, exportgraphics(gcf,fullfile(outdir,'mld_proxy.png'),'Resolution',300); end


fprintf('\n✅ All slab and sensitivity tests completed. Results saved to %s\n', outdir);


%% --- Supporting Function ---
function t50 = compute_tau50_rise(t, f)
    t50 = NaN;
    if isempty(t) || isempty(f) || all(~isfinite(f)) || max(f)<=0, return; end
    [fmax, imax] = max(f); 
    if imax<2, return; end
    target = 0.5*fmax;
    i = find(f(1:imax) >= target, 1, 'first');
    if isempty(i) || i==1, return; end
    t50 = t(i-1) + (target - f(i-1)) * (t(i) - t(i-1)) / (f(i) - f(i-1));
end