%% === SLAB MODEL: STRONGER VARIABLE ε + NONLINEAR BREAKUP (TEST 02) ===
%increases turbulence amplitude & nonlinearity
% --------------------------------------------------------------
clear; clc; close all;

outdir = fullfile(pwd, 'results_TEST3');
if ~exist(outdir, 'dir'), mkdir(outdir); end
save_figs = true;
colors = lines(8);
lw = 1.4; fs = 11;

%% === BASE CONFIGURATION ===
base_cfg = SimulationConfig();
base_cfg.t_final      = 30;          % duration [days]

% === Biological forcing (NPP): stronger ±80% ===
base_cfg.use_NPP      = true;
base_cfg.NPP_rate     = 5e-4;
base_cfg.NPP_profile  = 'sine';      % sinusoidal forcing (±80%)
% (will be plotted manually below with ±0.8 amplitude)

% === Turbulence forcing (ε): stronger + faster + phase-shifted ===
base_cfg.epsilon_profile = 'sine';
base_cfg.epsilon_mean    = 1e-6;
base_cfg.epsilon_amp     = 2.7e-6;   % ↑ amplitude (3× stronger than Test_02)
base_cfg.epsilon_period  = 5;        % ↓ shorter period (was 10 d)
base_cfg.epsilon_phase   = pi/2;     % 90° out of phase with NPP
base_cfg.epsilon         = base_cfg.epsilon_mean;  % compatibility

base_cfg.temp = 20;
base_cfg.salt = 35;

%% FIG 1 — Baseline (out-of-phase forcing)
fprintf('\n=== Running Test 03: Stronger ε forcing + nonlinear breakup ===\n');

cfg = base_cfg;
sim = CoagulationSimulation(cfg);
out = sim.run();
[t_flux, flux] = OutputGenerator.simpleFluxTimeSeries(out, DerivedGrid(cfg));

% === Plot time series ===
figure('Name','Test3_Baseline_StrongEps','Color','w');
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

tt = linspace(0, base_cfg.t_final, 200);
NPP_t = base_cfg.NPP_rate * (1 + 0.8*sin(2*pi*tt/base_cfg.t_final)); % ±80%
eps_t = base_cfg.epsilon_mean + base_cfg.epsilon_amp * ...
        sin(2*pi*tt/base_cfg.epsilon_period + base_cfg.epsilon_phase);

nexttile;
plot(tt, NPP_t, 'k', 'LineWidth', 1.5);
xlabel('Time (d)'); ylabel('NPP (model d⁻¹)');
title('NPP vs Time (±80 %)','FontWeight','bold'); grid on;

nexttile;
plot(tt, eps_t, 'r', 'LineWidth', 1.5);
xlabel('Time (d)'); ylabel('\epsilon (W kg⁻¹)');
title('\epsilon vs Time (90° out of phase)','FontWeight','bold');
grid on; ylim([min(eps_t)*0.8, max(eps_t)*1.2]);

nexttile;
flux_norm = flux / max(flux);
plot(t_flux, flux_norm, 'b', 'LineWidth', 1.5);
xlabel('Time (d)'); ylabel('Normalized Flux (a.u.)');
title('Flux vs Time (Normalized Response)','FontWeight','bold');
grid on; ylim([0 1.2]);

if save_figs
    exportgraphics(gcf, fullfile(outdir,'01_Test3_Baseline_StrongEps.png'),'Resolution',300);
end

fprintf('\n✅ Baseline figure saved (Test 03, stronger ε).\n');

%% FIG 2 — ε sweep (flux + τ₅₀)
fprintf('\n[1] Running ε sweep (strong forcing)…\n');
eps_list = [1e-8 1e-7 1e-6 1e-5];
tau50_eps = NaN(size(eps_list));
flux_curves = cell(numel(eps_list),1);
t_curves = cell(numel(eps_list),1);
Y_final = cell(numel(eps_list),1);
grid_all = cell(numel(eps_list),1);

for k = 1:numel(eps_list)
    cfgk = base_cfg;
    cfgk.epsilon_profile = 'constant';
    cfgk.epsilon = eps_list(k);
    simk = CoagulationSimulation(cfgk);
    outk = simk.run();
    gridk = DerivedGrid(cfgk);
    [tk, fk] = OutputGenerator.simpleFluxTimeSeries(outk, gridk);
    t_curves{k} = tk; flux_curves{k} = fk;
    tau50_eps(k) = compute_tau50_cumulative(tk, fk);
    Y_final{k} = outk.concentrations(end,:)';
    grid_all{k} = gridk;
end

figure('Name','Test3_eps_sweep','Color','w');
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on;
for k = 1:numel(eps_list)
    plot(t_curves{k}, flux_curves{k}, 'LineWidth', lw, 'Color', colors(k,:));
end
xlabel('Time (d)'); ylabel('Flux (a.u.)');
title('Flux vs Time — ε sweep (Test 03)','FontWeight','bold');
legend(arrayfun(@(e)sprintf('\\epsilon=%.0e',e),eps_list,'UniformOutput',false),'Location','best');
grid on;

nexttile;
mask = isfinite(tau50_eps);
if any(mask)
    semilogx(eps_list(mask), tau50_eps(mask),'o-k','LineWidth',lw,'MarkerFaceColor','w');
    grid on; xlabel('\epsilon (W kg⁻¹)'); ylabel('\tau_{50} (d)');
    title('\tau_{50} vs \epsilon (cumulative export)','FontWeight','bold');
    set(gca,'FontSize',fs);
else
    text(1e-7,0.5,'No valid τ₅₀','FontSize',fs);
end

if save_figs
    exportgraphics(gcf, fullfile(outdir,'02_Test3_eps_sweep.png'),'Resolution',300);
end

%% FIG 3 — Final PSDs for ε cases
figure('Name','Test3_PSD_vsEps','Color','w'); hold on;
nE = numel(eps_list); slopes = nan(1,nE);
for k = 1:nE
    if isempty(Y_final{k}), continue; end
    gridk = grid_all{k};
    diam = 2 * gridk.getFractalRadii() * 1e4; % cm→µm
    y = Y_final{k};
    valid = y>0 & isfinite(y);
    if nnz(valid)<3, continue; end
    D = diam(valid); N = y(valid);
    p = polyfit(log10(D), log10(N), 1);
    slopes(k) = -p(1);
    fit_line = 10.^(polyval(p, log10(D)));
    loglog(D,N,'LineWidth',1.5,'Color',colors(k,:));
    loglog(D,fit_line,'--','LineWidth',1.1,'Color',colors(k,:)*0.6);
end
xlabel('Particle diameter (µm)','FontSize',fs);
ylabel('Concentration (a.u.)','FontSize',fs);
title('Final PSD vs Turbulence (Test 03)','FontWeight','bold');
legend(arrayfun(@(e,b)sprintf('\\epsilon=%.0e (B=%.2f)',e,b),eps_list,slopes,'UniformOutput',false),...
       'Location','southwest','FontSize',fs);
set(gca,'XScale','log','YScale','log','XLim',[1e0 1e4],'YLim',[1e-10 1e-4],'LineWidth',1,'FontSize',fs);
grid on; box on;
if save_figs
    exportgraphics(gcf, fullfile(outdir,'03_Test3_PSD_vsEps.png'),'Resolution',300);
end
fprintf('\n✅ PSD figure with fitted slopes saved (Test 03).\n');

%% === Supporting function ===
function t50 = compute_tau50_cumulative(t,f)
    t50 = NaN;
    if isempty(t)||isempty(f)||all(~isfinite(f))||max(f)<=0, return; end
    t = t(:); f = f(:);
    f(~isfinite(f)|f<0)=0;
    Fcum = cumtrapz(t,f);
    Ftot = Fcum(end);
    if Ftot<=0, return; end
    target = 0.5*Ftot;
    i = find(Fcum>=target,1,'first');
    if isempty(i)||i==1
        t50 = t(1);
    else
        t50 = t(i-1)+(target-Fcum(i-1))*(t(i)-t(i-1))/(Fcum(i)-Fcum(i-1));
    end
end