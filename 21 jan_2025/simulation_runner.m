function simulation_runner
% SIMULATION_RUNNER (CONFIG-COMPATIBLE, STABLE VERSION)
%
% Runs:
%   1) Atlantic column (observed epsilon from epsilon_daily.mat)
%   2) Pulse turbulence experiment (baseline → pulse → baseline)  
%   3) EOF analysis (Atlantic run)
%
% Uses YOUR SimulationConfig fields:
%   - mu via cfg.growth_mode='pp' and cfg.growth
%   - disaggregation via cfg.enable_disagg + cfg.disagg_* fields
%
% Notes:
% - If disagg is ON, pulse runs can become stiff. This script:
%   (i) smooths pulses, (ii) sets safer ode options, (iii) asserts solver finished.

clc;

% -----------------------------
% USER TUNABLES
% -----------------------------
mu_tuned        = 0.12;       % [d^-1]

% Disaggregation switch (IMPORTANT: must be true/false, not "ture")
ENABLE_DISAGG   = true;       % <- set false for baseline

% Disagg tuning (start gentle)
DISAGG_KMAX_A   = 0.25;
DISAGG_BETA     = 0.10;
DISAGG_N_EXP    = 0.25;
DISAGG_PREDIST  = 0;

% Column settings
ZMAX = 200;   % [m]
DZ   = 10;    % [m]

% Output directory
outDir = 'fig_modelMain';
if ~exist(outDir,'dir'); mkdir(outDir); end

fprintf('======================================================\n');
fprintf('   SIMULATION RUNNER (Config-compatible)\n');
fprintf('   mu = %.3f d^-1 | disagg=%d\n', mu_tuned, ENABLE_DISAGG);
fprintf('======================================================\n');

%% ==============================================================
%  1) OBSERVED EPSILON (ATLANTIC COLUMN RUN)
% ===============================================================
[simA, cfgA, t_forcingA, eps_forcingA] = run_atlantic_column( ...
    mu_tuned, ZMAX, DZ, ENABLE_DISAGG, ...
    DISAGG_KMAX_A, DISAGG_BETA, DISAGG_N_EXP, DISAGG_PREDIST);

assert_solver_finished(simA, cfgA, 'Observed Atlantic run');

% NEW-2025-12-21: pass optional baseline window to plotting
baseline_windowA = [];  % let plotter auto-pick late-time stable window

figA = make_flux_plots(simA.result, simA.grid, cfgA, ...
    t_forcingA, eps_forcingA, ...
    sprintf('Observed (mu = %.2f, disagg=%d)', mu_tuned, ENABLE_DISAGG), ...
    baseline_windowA);

saveas(figA, fullfile(outDir, 'flux.png'));

%% ==============================================================
%  2) TURBULENCE PULSE EXPERIMENT (baseline → pulse → baseline)
% ===============================================================
[simB, cfgB, t_forcingB, eps_forcingB, metaB] = run_pulse_experiment( ...
    mu_tuned, ZMAX, DZ, ENABLE_DISAGG, ...
    DISAGG_KMAX_A, DISAGG_BETA, DISAGG_N_EXP, DISAGG_PREDIST);

assert_solver_finished(simB, cfgB, 'Pulse run');

% NEW-2025-12-21: define baseline window = pre-pulse
baseline_windowB = [0 metaB.t_pulse_start];

figB = make_flux_plots(simB.result, simB.grid, cfgB, ...
    t_forcingB, eps_forcingB, ...
    sprintf('Pulse (mu = %.2f, disagg=%d)', mu_tuned, ENABLE_DISAGG), ...
    baseline_windowB);

saveas(figB, fullfile(outDir, 'pulse_flux.png'));

%% ==============================================================
%  3) EOF ANALYSIS (ATLANTIC RUN ONLY)
% ===============================================================
tA    = simA.result.time(:);
YA    = simA.result.concentrations;
gridA = simA.grid;
D_um  = 2 * gridA.getFractalRadii() * 1e4;  % diameter in microns

[figC_maps, figC_3d] = analyze_model_eofs_save( ...
    'Atlantic', tA, YA, D_um, cfgA, ...
    t_forcingA, eps_forcingA, gridA);

saveas(figC_maps, fullfile(outDir, 'EOF_maps_PCs.png'));
saveas(figC_3d,   fullfile(outDir, 'EOF_3D_slices.png'));

fprintf('\nDONE. Figures saved in: %s\n', outDir);

end

%% =================================================================
% ASSERT: SOLVER DID NOT STOP EARLY
%% =================================================================
function assert_solver_finished(simObj, cfg, label)
t = simObj.result.time(:);
t_end = t(end);
t_goal = cfg.t_final;

% allow tiny tolerance
if t_end < (t_goal - 1e-6)
    error('%s: solver stopped early at t=%.4f (expected %.4f). Reduce disagg strength or tighten ODE settings.', ...
        label, t_end, t_goal);
end
end

%% =================================================================
% 1) COLUMN RUN WITH OBSERVED EPSILON
%% =================================================================
function [sim, cfg, t_forcing, eps_forcing] = run_atlantic_column( ...
    mu_val, z_max, dz, enable_disagg, kmax_a, beta, n_exp, predis)

if ~isfile('epsilon_daily.mat')
    error('Missing epsilon_daily.mat in current folder.');
end

S = load('epsilon_daily.mat');
[t_days, eps_raw] = grab_data_simple(S);

mask        = t_days <= 60;
t_forcing   = t_days(mask);
eps_forcing = eps_raw(mask);

t_forcing   = t_forcing(:) - t_forcing(1);
eps_forcing = eps_forcing(:);

fprintf('Running Atlantic column for %.1f days (N=%d forcing steps)...\n', ...
    range(t_forcing), numel(t_forcing));

cfg = SimulationConfig();
cfg.debug_rhs_units = false;   % NEW-2025-12-21: keep OFF for normal runs
cfg.use_column  = true;
cfg.z_max       = z_max;
cfg.dz          = dz;

cfg.growth_mode = 'pp';
cfg.growth      = mu_val;

cfg.t_init  = t_forcing(1);
cfg.t_final = t_forcing(end);
if numel(t_forcing) > 1
    cfg.delta_t = min(diff(t_forcing));
end

% Disaggregation config
cfg.enable_disagg   = logical(enable_disagg);
cfg.disagg_apply_in = 'rhs';
cfg.disagg_operator = 'applyWithScaling';

cfg.eps_ref               = 1e-6;
cfg.disagg_kmax_a         = kmax_a;
cfg.disagg_beta           = beta;
cfg.disagg_n_exp          = n_exp;
cfg.disagg_redistribute_p = predis;

cfg.epsilon_time   = t_forcing;
cfg.epsilon_series = eps_forcing;

% ODE stability upgrades (especially when disagg is ON)
cfg.clip_negative  = true;
if cfg.enable_disagg
    cfg.use_nonnegative = false;  % NonNegative can make stiff events harder
    cfg.ode_options = odeset('RelTol',1e-6,'AbsTol',1e-12,'MaxStep',0.05,'InitialStep',1e-4);
else
    cfg.ode_options = odeset('RelTol',1e-6,'AbsTol',1e-12);
end

sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_forcing, eps_forcing);
sim.run('tspan', t_forcing);

end

%% =================================================================
% 2) TURBULENCE PULSE EXPERIMENT (baseline → pulse → baseline)
%% =================================================================
function [sim, cfg, t_forcing, eps_forcing, meta] = run_pulse_experiment( ...
    mu_val, z_max, dz, enable_disagg, kmax_a, beta, n_exp, predis)

% ==============================================================

%   baseline for 20 days
%   pulse for 1 day (or change to 5)
%   recovery for 10 days
% ==============================================================
meta = struct();
meta.T_baseline   = 20;   % days
meta.T_pulse      = 1;    % days  (try 5 later)
meta.T_recovery   = 10;   % days
meta.dt           = 0.05; % days

meta.t_pulse_start = meta.T_baseline;
meta.t_pulse_end   = meta.T_baseline + meta.T_pulse;

T_total  = meta.T_baseline + meta.T_pulse + meta.T_recovery;
dt       = meta.dt;
t_forcing = (0:dt:T_total).';

eps_base  = 1e-9;

% NEW-2025-12-21: choose one pulse strength now (you can sweep later)
eps_peak  = 5e-6;

% Smooth pulse using tanh ramps (stable for ode15s)
ramp = 0.05; % days
t1 = meta.t_pulse_start;
t2 = meta.t_pulse_end;

up     = 0.5*(1 + tanh((t_forcing - t1)/ramp));
down   = 0.5*(1 + tanh((t2 - t_forcing)/ramp));
window = up .* down;

eps_forcing = eps_base * ones(size(t_forcing));
eps_forcing = eps_forcing + (eps_peak - eps_base) * window;

fprintf('Running Pulse experiment for %.1f days (N=%d forcing steps)...\n', ...
    range(t_forcing), numel(t_forcing));

cfg = SimulationConfig();
cfg.debug_rhs_units = false;   % NEW-2025-12-21: keep OFF for normal runs
cfg.use_column  = true;
cfg.z_max       = z_max;
cfg.dz          = dz;

cfg.growth_mode = 'pp';
cfg.growth      = mu_val;

cfg.t_init   = t_forcing(1);
cfg.t_final  = t_forcing(end);
cfg.delta_t  = dt;

cfg.enable_disagg   = logical(enable_disagg);
cfg.disagg_apply_in = 'rhs';
cfg.disagg_operator = 'applyWithScaling';

cfg.eps_ref               = 1e-6;
cfg.disagg_kmax_a         = kmax_a;
cfg.disagg_beta           = beta;
cfg.disagg_n_exp          = n_exp;
cfg.disagg_redistribute_p = predis;

cfg.epsilon_time   = t_forcing;
cfg.epsilon_series = eps_forcing;

cfg.clip_negative  = true;
if cfg.enable_disagg
    cfg.use_nonnegative = false;
    cfg.ode_options = odeset('RelTol',1e-6,'AbsTol',1e-12,'MaxStep',0.02,'InitialStep',1e-5);
else
    cfg.ode_options = odeset('RelTol',1e-6,'AbsTol',1e-12);
end

sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_forcing, eps_forcing);
sim.run('tspan', t_forcing);

end

%% =================================================================
% FLUX + SIZE-CLASS PLOTS
%% =================================================================
function fig = make_flux_plots(result, grid_obj, cfg, ...
    t_forcing, eps_forcing, title_prefix, baseline_window)

t = result.time(:);
Y = result.concentrations;

Ns = cfg.n_sections;
Nz = cfg.getNumLayers();

eps_t = interp1(t_forcing(:), eps_forcing(:), t, 'previous', 'extrap');
eps_t = max(eps_t, 1e-12);

r_cm    = grid_obj.getFractalRadii();
r_v     = grid_obj.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_obj.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;
V_m3    = (4/3) * pi * r_v.^3 * 1e-6;
D_um    = 2 * r_cm * 1e4;

export_depth_idx = Nz;
c1 = (export_depth_idx-1)*Ns + 1;
c2 = export_depth_idx*Ns;
Q_export = Y(:, c1:c2);

F_total = sum(Q_export .* V_m3.' .* W_m_d.', 2);

% -------------------------------------------------------------
% OLD (kept):
% idx0 = find(eps_t < 1e-10, min(10,numel(t)));
% if isempty(idx0), idx0 = 1:min(10,numel(t)); end
% F_ss = mean(F_total(idx0));
% if ~isfinite(F_ss) || F_ss == 0, F_ss = 1; end
% F_relative = F_total / F_ss;
% -------------------------------------------------------------

% ==========================================================
% NEW-2025-12-21: robust baseline for relative flux
% If baseline_window is given [t0 t1], compute F_ss from that window.
% Else fallback to late-time median.
% ==========================================================
idx_ref = [];

if nargin >= 7 && ~isempty(baseline_window) && numel(baseline_window) == 2
    idx_ref = find(t >= baseline_window(1) & t <= baseline_window(2));
    % use last part of baseline to avoid startup transient
    if numel(idx_ref) > 50
        idx_ref = idx_ref(end-49:end);
    end
end

if isempty(idx_ref)
    % fallback: last 50 points
    idx_ref = max(1, numel(t)-49):numel(t);
end

F_ss = median(F_total(idx_ref));          % robust to spikes / transients
if ~isfinite(F_ss) || abs(F_ss) < 1e-30
    F_ss = 1e-30;
end

F_relative = F_total / F_ss;

% OLD:
% D_small_max  = 200;
% D_medium_max = 1000;

% NEW-2026-01-06 
% small < 500 µm
% medium 500–2000 µm
% large > 2000 µm
D_small_max  = 500;
D_medium_max = 2000;

idx_small  = D_um < D_small_max;
idx_medium = (D_um >= D_small_max) & (D_um < D_medium_max);
idx_large  = D_um >= D_medium_max;

F_small  = sum(Q_export(:, idx_small)  .* V_m3(idx_small).'  .* W_m_d(idx_small).', 2);
F_medium = sum(Q_export(:, idx_medium) .* V_m3(idx_medium).' .* W_m_d(idx_medium).', 2);
F_large  = sum(Q_export(:, idx_large)  .* V_m3(idx_large).'  .* W_m_d(idx_large).', 2);

F_total_safe = max(F_total, 1e-30);
frac_small   = F_small  ./ F_total_safe;
frac_medium  = F_medium ./ F_total_safe;
frac_large   = F_large  ./ F_total_safe;

fig = figure('Name',[title_prefix ' Flux Analysis'], ...
    'Color','w','Position',[100 100 800 900]);

subplot(3,1,1);
semilogy(t, eps_t, 'k-', 'LineWidth',1.2); hold on;
plot(t([1 end]), [1e-6 1e-6], 'r--', 'LineWidth',1);
ylabel('\epsilon (W kg^{-1})');
title(title_prefix);
axis tight; grid on;

subplot(3,1,2);
plot(t, F_relative, 'b-', 'LineWidth',2); hold on;
yline(1,'k--','LineWidth',1);
ylabel('Relative Flux (F/F_{ss})');
title('Relative Export Flux');
axis tight; grid on;

subplot(3,1,3);
plot(t, frac_large,  'r-','LineWidth',2.0,'DisplayName','Large'); hold on;
plot(t, frac_medium,'g-','LineWidth',2.0,'DisplayName','Medium');
plot(t, frac_small, 'b-','LineWidth',2.0,'DisplayName','Small');
ylim([0 1]);
ylabel('Flux Fraction');
xlabel('Time (days)');
legend('Location','northeast');
grid on; box on;

end

%% =================================================================
% EOF ANALYSIS (upgraded to avoid PCA TSQUARED warning)
% (UNCHANGED from your version below)
%% =================================================================
function [fig_maps, fig_3d] = analyze_model_eofs_save( ...
    title_prefix, t, Y, D_um, cfg, t_forcing, eps_forcing, grid_obj)

t = t(:);
Nt = numel(t);
Ns = cfg.n_sections;
Nz = cfg.getNumLayers();
depths = cfg.getZ(); depths = depths(:);

eps_t = interp1(t_forcing(:), eps_forcing(:), t, 'previous', 'extrap');
eps_t = max(eps_t, 1e-13);

Y_3d = zeros(Nt, Nz, Ns);
for k = 1:Nz
    c1 = (k-1)*Ns + 1;
    c2 = k*Ns;
    Y_3d(:,k,:) = Y(:, c1:c2);
end

Y_flat = reshape(Y_3d, Nt, Nz*Ns);
Y_log  = log10(max(Y_flat, 1e-20));
Y_mean = mean(Y_log, 1);
Y_anom = Y_log - Y_mean;

col_std = std(Y_anom, 0, 1);
keep = isfinite(col_std) & (col_std > 1e-12);
if ~any(keep)
    keep = true(1, size(Y_anom,2));
end
Y_anom2 = Y_anom(:, keep);

coeffs2 = [];
scores2 = [];
explained = [];
try
    [coeffs2, scores2, latent] = pca(Y_anom2, ...
        'Algorithm','svd','Economy',true,'Centered',false);
    explained = 100 * latent(:)' / sum(latent);
catch
    [U,S,V] = svd(Y_anom2, 'econ');
    coeffs2 = V;
    scores2 = U*S;
    svals = diag(S);
    lam = svals.^2;
    explained = 100 * lam(:)' / sum(lam);
end

coeffs = zeros(size(Y_anom,2), size(coeffs2,2));
scores = scores2;
coeffs(keep, :) = coeffs2;

D_mm = D_um(:) / 1000;
logD = log10(D_mm);

mode1 = reshape(coeffs(:,1), Ns, Nz).';
mode2 = reshape(coeffs(:,2), Ns, Nz).';

if mean(mode1(1,end-3:end)) < 0
    mode1       = -mode1;
    coeffs(:,1) = -coeffs(:,1);
    scores(:,1) = -scores(:,1);
end
if mean(mode2(1,end-3:end)) < 0
    mode2       = -mode2;
    coeffs(:,2) = -coeffs(:,2);
    scores(:,2) = -scores(:,2);
end

fig_maps = figure('Name',[title_prefix ' EOF Analysis'], ...
    'Color','w','Position',[50 50 1200 850]);

subplot(2,2,1);
semilogy(t, eps_t, 'k-','LineWidth',1.2);
ylabel('\epsilon (W kg^{-1})');
xlabel('Time (days)');
title('Turbulence Forcing');
axis tight; grid on;

subplot(2,2,2);
pc1 = scores(:,1); pc2 = scores(:,2);
pc1_n = pc1 ./ max(abs(pc1));
pc2_n = pc2 ./ max(abs(pc2));

export_depth_idx = Nz;
c1 = (export_depth_idx-1)*Ns + 1;
c2 = export_depth_idx*Ns;
Q_export = Y(:, c1:c2);

r_cm    = grid_obj.getFractalRadii();
r_v     = grid_obj.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_obj.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;
V_m3    = (4/3) * pi * r_v.^3 * 1e-6;

F_total = sum(Q_export .* V_m3.' .* W_m_d.', 2);
F_rel   = F_total / max(F_total);

plot(t, pc1_n, 'k-','LineWidth',1.2,'DisplayName','PC1'); hold on;
plot(t, pc2_n, 'r-','LineWidth',1.2,'DisplayName','PC2');
plot(t, F_rel, 'k--','LineWidth',1.0,'DisplayName','Flux (norm)');
yline(0,'k:');
ylabel('Normalized units');
xlabel('Time (days)');
title(sprintf('PC1 & PC2 (%.1f%% / %.1f%% var)', explained(1), explained(2)));
legend('Location','best');
axis tight; grid on;

subplot(2,2,3);
imagesc(logD, depths, mode1);
set(gca,'YDir','reverse');
colormap(gca, turbo);
cb = colorbar; cb.Label.String = 'EOF1 loading';
title(sprintf('Mode 1 (%.1f%%%% var)', explained(1)));
xlabel('D (mm), log_{10}');
ylabel('Depth (m)');
set(gca,'XTick',log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});

subplot(2,2,4);
imagesc(logD, depths, mode2);
set(gca,'YDir','reverse');
colormap(gca, turbo);
cb = colorbar; cb.Label.String = 'EOF2 loading';
title(sprintf('Mode 2 (%.1f%%%% var)', explained(2)));
xlabel('D (mm), log_{10}');
ylabel('Depth (m)');
set(gca,'XTick',log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});

PSD_full_3d = reshape(Y_log, Nt, Nz, Ns);

Y1_anom = scores(:,1) * coeffs(:,1).';
Y1_log  = Y1_anom + Y_mean;
Y1_3d   = reshape(Y1_log, Nt, Nz, Ns);

Y2_anom = scores(:,2) * coeffs(:,2).';
Y2_log  = Y2_anom + Y_mean;
Y2_3d   = reshape(Y2_log, Nt, Nz, Ns);

vals = sort(Y_log(:));
vals = vals(isfinite(vals));
if isempty(vals)
    clim_global = [-20 -5];
else
    n  = numel(vals);
    lo = vals(max(1, round(0.02*n)));
    hi = vals(round(0.98*n));
    if hi <= lo, hi = lo + 1; end
    clim_global = [lo hi];
end

fig_3d = figure('Name',[title_prefix ' 3D slices'], ...
    'Color','w','Position',[50 50 1100 950]);

ax1 = subplot(3,1,1);
plot_psd_slice_subplot(ax1, t, depths, D_mm, PSD_full_3d, ...
    'Full Model (log_{10})', clim_global);

ax2 = subplot(3,1,2);
plot_psd_slice_subplot(ax2, t, depths, D_mm, Y1_3d, ...
    'Mode 1 reconstruction', clim_global);

ax3 = subplot(3,1,3);
plot_psd_slice_subplot(ax3, t, depths, D_mm, Y2_3d, ...
    'Mode 2 reconstruction', clim_global);

end

function plot_psd_slice_subplot(ax_handle, t_days, depths, D_mm, PSD_3d, title_str, c_limits)
axes(ax_handle);
hold on;

t_days = t_days(:);
depths = depths(:);
logD   = log10(D_mm(:));

V = permute(PSD_3d, [3 1 2]);   % [Ns x Nt x Nz]
[X,Y,Z] = meshgrid(t_days, logD, depths);

t_min = ceil(min(t_days));
t_max = floor(max(t_days));
xs    = t_min:1:t_max;

h = slice(X,Y,Z,V,fliplr(xs),[],[]);
set(h,'EdgeColor','none','FaceAlpha',0.95);

set(gca,'ZDir','reverse');
colormap(gca, turbo);

if numel(c_limits)==2 && c_limits(2) > c_limits(1)
    caxis(c_limits);
end

cb = colorbar('Location','eastoutside');
cb.Label.String = 'log_{10} Part Vol';

xlabel('Time (days)');
ylabel('D (mm), log_{10}');
zlabel('Depth (m)');
title(title_str);

set(gca,'YTick',log10([0.1 1 10]), 'YTickLabel',{'0.1','1','10'});
pbaspect([35 1 3]);
view([-60 5]);
axis tight; box on; grid on;
end

function [t, e] = grab_data_simple(S)
names = fieldnames(S);
if numel(names)==1 && isstruct(S.(names{1}))
    S     = S.(names{1});
    names = fieldnames(S);
end

t = [];
e = [];
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

t = t(:);
e = e(:);

t = t - t(1);
n = min(numel(t), numel(e));
t = t(1:n);
e = e(1:n);

e = max(e, 1e-12);
end