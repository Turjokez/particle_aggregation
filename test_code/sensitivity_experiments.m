%% === #example_slab_sensitivity_experiments.m ===

clear; clc; close all;

save_figs = true;
outdir = fullfile(pwd,'results');
if save_figs && ~exist(outdir,'dir'), mkdir(outdir); end
colors = lines(6);

%% τ50 (time to half-peak on rising limb)
tau50_rise = @(t,f) compute_tau50_rise(t,f);

%% Base configuration 
cfg = SimulationConfig();
cfg.t_final     = 30;
cfg.use_NPP     = true;
cfg.NPP_rate    = 5e-4;
cfg.NPP_profile = 'sine';
cfg.epsilon     = 1e-6;
cfg.temp = 20; cfg.salt = 35;

lw = 1.0; fs = 11;    % thin lines, readable text

%% (1) BASELINE TEST 
fprintf('\n[1] Running baseline slab test...\n');
tt = linspace(0, cfg.t_final, 100);
NPP_t = cfg.NPP_rate * (1 + 0.5*sin(2*pi*tt/cfg.t_final));
eps_t = cfg.epsilon * ones(size(tt));

sim  = CoagulationSimulation(cfg);
out  = sim.run();
grid = DerivedGrid(cfg);
[~, flux] = OutputGenerator.simpleFluxTimeSeries(out, grid);

figure('Name','Core_Baseline','Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');

nexttile; plot(tt, NPP_t, 'b','LineWidth',lw);
xlabel('Time (d)','FontSize',fs); ylabel('NPP (model d^{-1})','FontSize',fs);
title('NPP vs Time','FontWeight','bold'); grid on;

nexttile; semilogy(tt, eps_t, 'r','LineWidth',lw);
xlabel('Time (d)','FontSize',fs); ylabel('\epsilon (W kg^{-1})','FontSize',fs);
title('\epsilon vs Time','FontWeight','bold'); grid on; ylim([1e-8 1e-5]);

nexttile; plot(out.time, flux, 'k','LineWidth',lw);
xlabel('Time (d)','FontSize',fs); ylabel('Flux (a.u.)','FontSize',fs);
title('Flux vs Time','FontWeight','bold'); grid on;

if save_figs, exportgraphics(gcf, fullfile(outdir,'Core_Baseline.png'),'Resolution',300); end


%% (2) EPSILON SWEEP 
fprintf('\n[2] Running ε sweep...\n');
eps_list = [1e-8 1e-7 1e-6 1e-5];
tau50_eps = NaN(size(eps_list));
flux_curves = cell(numel(eps_list),1);
t_curves    = cell(numel(eps_list),1);

for k = 1:numel(eps_list)
    cfgk = cfg;
    cfgk.epsilon = eps_list(k);
    simk  = CoagulationSimulation(cfgk);
    outk  = simk.run();
    gridk = DerivedGrid(cfgk);
    [tk, fk] = OutputGenerator.simpleFluxTimeSeries(outk, gridk);
    t_curves{k} = tk;
    flux_curves{k} = fk;
    tau50_eps(k) = tau50_rise(tk,fk);
end

figure('Name','eps_sweep','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% --- Left: Flux vs Time ---
nexttile; hold on;
for k=1:numel(eps_list)
    plot(t_curves{k}, flux_curves{k}, 'LineWidth',lw, 'Color',colors(k,:));
end
hold off; grid on;
xlabel('Time (d)','FontSize',fs); ylabel('Flux (a.u.)','FontSize',fs);
title('Flux vs Time — \epsilon sweep','FontWeight','bold');
legend(arrayfun(@(e) sprintf('\\epsilon = %.0e', e), eps_list,'UniformOutput',false), 'Location','best');

% --- Right: τ50 vs ε ---
nexttile; mask = isfinite(tau50_eps) & tau50_eps > 0;
if any(mask)
    semilogx(eps_list(mask), tau50_eps(mask), 'o-','LineWidth',lw*1.2,...
        'Color','k','MarkerFaceColor','w');
else
    text(1e-7,0.5,'No valid \tau_{50} values','FontSize',fs);
end
grid on; xlabel('\epsilon (W kg^{-1})','FontSize',fs);
ylabel('\tau_{50} (d)','FontSize',fs);
title('\tau_{50} vs \epsilon','FontWeight','bold');
if save_figs, exportgraphics(gcf, fullfile(outdir,'eps_sweep.png'),'Resolution',300); end


%% (3) NPP RESPONSE
fprintf('\n[3] Running NPP sensitivity & step test...\n');
npp_levels = [1e-4 5e-4 1e-3];
tN = cell(numel(npp_levels),1); fN = cell(numel(npp_levels),1);

for k=1:numel(npp_levels)
    cfgk = cfg;
    cfgk.NPP_rate = npp_levels(k);
    cfgk.NPP_profile = 'constant';
    outk = CoagulationSimulation(cfgk).run();
    [tN{k}, fN{k}] = OutputGenerator.simpleFluxTimeSeries(outk, DerivedGrid(cfgk));
end

% Step change test
cfg_stepA = cfg; cfg_stepA.t_final = 15; cfg_stepA.NPP_profile='constant'; cfg_stepA.NPP_rate = 5e-4;
outA = CoagulationSimulation(cfg_stepA).run();
[tA, fA] = OutputGenerator.simpleFluxTimeSeries(outA, DerivedGrid(cfg_stepA));

cfg_stepB = cfg; cfg_stepB.t_final = 15; cfg_stepB.NPP_profile='constant'; cfg_stepB.NPP_rate = 1e-3;
outB = CoagulationSimulation(cfg_stepB).run();
[tB, fB] = OutputGenerator.simpleFluxTimeSeries(outB, DerivedGrid(cfg_stepB));

t_step = [tA; tB+15]; f_step = [fA; fB];

figure('Name','npp_response','Color','w');
tiledlayout(1,2,'Padding','compact','TileSpacing','compact');

% constant NPP levels
nexttile; hold on;
for k=1:numel(npp_levels)
    plot(tN{k}, fN{k}, 'LineWidth',lw, 'Color',colors(k,:));
end
hold off; grid on;
xlabel('Time (d)','FontSize',fs); ylabel('Flux (a.u.)','FontSize',fs);
title('Flux vs Time — NPP levels','FontWeight','bold');
legend(arrayfun(@(x) sprintf('NPP = %.1e',x), npp_levels,'UniformOutput',false),'Location','best');

% step test
nexttile;
plot(t_step, f_step, 'k','LineWidth',lw); grid on;
xlabel('Time (d)','FontSize',fs); ylabel('Flux (a.u.)','FontSize',fs);
title('Flux response to NPP step (day 15)','FontWeight','bold');
if save_figs, exportgraphics(gcf, fullfile(outdir,'npp_response.png'),'Resolution',300); end


%% (4) MLD PROXY TEST 
fprintf('\n[4] Running MLD proxy (10–50 m)...\n');
mld_list = [10 20 30 40 50]; ref_mld = 30;
tM = cell(numel(mld_list),1); fM = cell(numel(mld_list),1);

for k=1:numel(mld_list)
    cfgk = cfg;
    scale = ref_mld / mld_list(k);
    cfgk.NPP_rate = cfg.NPP_rate * scale;
    cfgk.NPP_profile = 'constant';
    outk = CoagulationSimulation(cfgk).run();
    [tM{k}, fM{k}] = OutputGenerator.simpleFluxTimeSeries(outk, DerivedGrid(cfgk));
end

figure('Name','mld_proxy','Color','w');
hold on;
for k=1:numel(mld_list)
    plot(tM{k}, fM{k}, 'LineWidth',lw, 'Color',colors(k,:));
end
hold off; grid on;
xlabel('Time (d)','FontSize',fs); ylabel('Flux (a.u.)','FontSize',fs);
title('Flux vs Time — mixed-layer proxy (10–50 m)','FontWeight','bold');
legend(arrayfun(@(z) sprintf('MLD %2.0f m',z), mld_list,'UniformOutput',false),'Location','best');
if save_figs, exportgraphics(gcf, fullfile(outdir,'mld_proxy.png'),'Resolution',300); end


%% Supporting function 
function t50 = compute_tau50_rise(t, f)
    % Compute time to 50% of max flux (rising limb only)
    t50 = NaN;
    if isempty(t) || isempty(f) || all(~isfinite(f)) || max(f) <= 0
        return
    end
    [fmax, imax] = max(f);
    if imax < 2, return; end
    target = 0.5 * fmax;
    i = find(f(1:imax) >= target, 1, 'first');
    if isempty(i) || i == 1, return; end
    % Linear interpolate between f(i-1) and f(i)
    t50 = t(i-1) + (target - f(i-1)) * (t(i) - t(i-1)) / (f(i) - f(i-1));
end