function simulation_runner
% SIMULATION_RUNNER (FULL UPDATED VERSION)
% - Explicit PP (enable_pp) + biovolume state bookkeeping
% - Flux plots: bigger + clean export + NO toolbar warning
% - EOF maps: diverging colormap with white at 0 (symmetric caxis)
% - size spectra plots (log N vs log D) at chosen depths/times

clc;
% Turn off axes toolbar globally (prevents exportgraphics warning)
set(groot,'defaultAxesToolbarVisible','off');
% -----------------------------
% USER TUNABLES
% -----------------------------
PP_AREAL      = 0.1;     % [cm^3/cm^2/day] because state_is_biovolume=true
ENABLE_DISAGG = false;   % keep OFF until baseline is trusted

% Disaggregation tuning (gentle)
DISAGG_KMAX_A  = 0.25;
DISAGG_BETA    = 0.10;
DISAGG_N_EXP   = 0.25;
DISAGG_PREDIST = 0;

% Column settings
ZMAX = 200;   % [m]
DZ   = 10;    % [m]

% Output directory
outDir = 'fig_modelMain';
if ~exist(outDir,'dir'); mkdir(outDir); end

fprintf('======================================================\n');
fprintf('   SIMULATION RUNNER (FULL UPDATED)\n');
fprintf('   PP_AREAL = %.4g | disagg=%d\n', PP_AREAL, ENABLE_DISAGG);
fprintf('======================================================\n');

%% ==============================================================
%  1) OBSERVED EPSILON (ATLANTIC COLUMN RUN)
% ===============================================================
[simA, cfgA, t_forcingA, eps_forcingA] = run_atlantic_column( ...
    PP_AREAL, ZMAX, DZ, ENABLE_DISAGG, ...
    DISAGG_KMAX_A, DISAGG_BETA, DISAGG_N_EXP, DISAGG_PREDIST);

assert_solver_finished(simA, cfgA, 'Observed Atlantic run');

baseline_windowA = []; % auto: late-time median window

figA = make_flux_plots(simA.result, simA.grid, cfgA, ...
    t_forcingA, eps_forcingA, ...
    sprintf('Observed (PP=%.3g, disagg=%d)', PP_AREAL, ENABLE_DISAGG), ...
    baseline_windowA);

exportgraphics(figA, fullfile(outDir,'flux.png'), 'Resolution', 200);

%% ==============================================================
%  2) TURBULENCE PULSE EXPERIMENT (baseline → pulse → baseline)
% ===============================================================
[simB, cfgB, t_forcingB, eps_forcingB, metaB] = run_pulse_experiment( ...
    PP_AREAL, ZMAX, DZ, ENABLE_DISAGG, ...
    DISAGG_KMAX_A, DISAGG_BETA, DISAGG_N_EXP, DISAGG_PREDIST);

assert_solver_finished(simB, cfgB, 'Pulse run');

baseline_windowB = [0 metaB.t_pulse_start];

figB = make_flux_plots(simB.result, simB.grid, cfgB, ...
    t_forcingB, eps_forcingB, ...
    sprintf('Pulse (PP=%.3g, disagg=%d)', PP_AREAL, ENABLE_DISAGG), ...
    baseline_windowB);

exportgraphics(figB, fullfile(outDir,'pulse_flux.png'), 'Resolution', 200);

%% ==============================================================
%  2b) size spectra plots (Observed + Pulse)
% ===============================================================
% Pick depths close to his suggestion (adjust if your z_max changes)
depths_m = [10 50 150];

% Times: early + around turbulence peak region
times_obs = [0 5 10 19 20 21 22 23 24 25];
times_pul = [0 5 10 18 19 20 20.5 21 22 25];

figS1 = plot_size_spectra_selected(simA.result, cfgA, ...
    'Observed: size spectra (lines = times)', depths_m, times_obs);
exportgraphics(figS1, fullfile(outDir,'Observed_sizeSpectra.png'), 'Resolution', 200);

figS2 = plot_size_spectra_selected(simB.result, cfgB, ...
    'Pulse: size spectra (lines = times)', depths_m, times_pul);
exportgraphics(figS2, fullfile(outDir,'Pulse_sizeSpectra.png'), 'Resolution', 200);

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

exportgraphics(figC_maps, fullfile(outDir,'EOF_maps_PCs.png'), 'Resolution', 200);
exportgraphics(figC_3d,   fullfile(outDir,'EOF_3D_slices.png'), 'Resolution', 200);

fprintf('\nDONE. Figures saved in: %s\n', outDir);

end

%% =================================================================
% ASSERT: SOLVER DID NOT STOP EARLY
%% =================================================================
function assert_solver_finished(simObj, cfg, label)
t = simObj.result.time(:);
t_end  = t(end);
t_goal = cfg.t_final;
if t_end < (t_goal - 1e-6)
    error('%s: solver stopped early at t=%.4f (expected %.4f).', label, t_end, t_goal);
end
end

%% =================================================================
% 1) COLUMN RUN WITH OBSERVED EPSILON
%% =================================================================
function [sim, cfg, t_forcing, eps_forcing] = run_atlantic_column( ...
    PP_AREAL, z_max, dz, enable_disagg, kmax_a, beta, n_exp, predis)

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
cfg.debug_rhs_units = false;

cfg.use_column = true;
cfg.z_max      = z_max;
cfg.dz         = dz;

% ---- PHASE-1 bookkeeping ----
cfg.state_is_biovolume = true;
cfg.export_weight      = "ones";
cfg.enable_sinking     = true;

% ---- explicit PP injection ----
cfg.growth_mode = 'shift';
cfg.growth      = 0;

cfg.enable_pp = true;
cfg.pp_rate   = PP_AREAL;   % cm^3/cm^2/day (biovolume state)
cfg.pp_bin    = 1;
cfg.pp_layer  = 1;

% Time
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

% ODE options
cfg.clip_negative = true;
if cfg.enable_disagg
    cfg.use_nonnegative = false;
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
    PP_AREAL, z_max, dz, enable_disagg, kmax_a, beta, n_exp, predis)

meta = struct();
meta.T_baseline = 20;
meta.T_pulse    = 1;
meta.T_recovery = 10;
meta.dt         = 0.05;

meta.t_pulse_start = meta.T_baseline;
meta.t_pulse_end   = meta.T_baseline + meta.T_pulse;

T_total   = meta.T_baseline + meta.T_pulse + meta.T_recovery;
dt        = meta.dt;
t_forcing = (0:dt:T_total).';

eps_base = 1e-9;
eps_peak = 5e-6;

% Smooth pulse (tanh ramps)
ramp = 0.05;
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
cfg.debug_rhs_units = false;

cfg.use_column = true;
cfg.z_max      = z_max;
cfg.dz         = dz;

% ---- PHASE-1 bookkeeping ----
cfg.state_is_biovolume = true;
cfg.export_weight      = "ones";
cfg.enable_sinking     = true;

% ---- explicit PP injection ----
cfg.growth_mode = 'shift';
cfg.growth      = 0;

cfg.enable_pp = true;
cfg.pp_rate   = PP_AREAL;
cfg.pp_bin    = 1;
cfg.pp_layer  = 1;

% Time
cfg.t_init  = t_forcing(1);
cfg.t_final = t_forcing(end);
cfg.delta_t = dt;

% Disaggregation
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

cfg.clip_negative = true;
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
% FLUX + SIZE-CLASS PLOTS (BIGGER + no toolbar + biovolume-safe)
%% =================================================================
function fig = make_flux_plots(result, grid_obj, cfg, ...
    t_forcing, eps_forcing, title_prefix, baseline_window)

t = result.time(:);
Y = result.concentrations;

Ns = cfg.n_sections;
Nz = cfg.getNumLayers();

eps_t = interp1(t_forcing(:), eps_forcing(:), t, 'previous', 'extrap');
eps_t = max(eps_t, 1e-14);

r_cm    = grid_obj.getFractalRadii();
r_v     = grid_obj.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_obj.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;

V_m3 = (4/3) * pi * r_v.^3 * 1e-6;  % used only if not biovolume-state
D_um = 2 * r_cm * 1e4;

% Bottom layer export
export_depth_idx = Nz;
c1 = (export_depth_idx-1)*Ns + 1;
c2 = export_depth_idx*Ns;
Q_export = Y(:, c1:c2);

% Export flux
if isprop(cfg,'state_is_biovolume') && logical(cfg.state_is_biovolume)
    F_total = sum(Q_export .* W_m_d.', 2);          % biovolume flux
else
    F_total = sum(Q_export .* V_m3.' .* W_m_d.', 2);
end

% Robust baseline for relative flux
idx_ref = [];
if nargin >= 7 && ~isempty(baseline_window) && numel(baseline_window) == 2
    idx_ref = find(t >= baseline_window(1) & t <= baseline_window(2));
    if numel(idx_ref) > 50, idx_ref = idx_ref(end-49:end); end
end
if isempty(idx_ref)
    idx_ref = max(1, numel(t)-49):numel(t);
end

F_ss = median(F_total(idx_ref));
if ~isfinite(F_ss) || abs(F_ss) < 1e-30
    F_ss = 1e-30;
end
F_relative = F_total / F_ss;

% bins
D_small_max  = 500;
D_medium_max = 2000;

idx_small  = D_um < D_small_max;
idx_medium = (D_um >= D_small_max) & (D_um < D_medium_max);
idx_large  = D_um >= D_medium_max;

if isprop(cfg,'state_is_biovolume') && logical(cfg.state_is_biovolume)
    F_small  = sum(Q_export(:, idx_small)  .* W_m_d(idx_small).', 2);
    F_medium = sum(Q_export(:, idx_medium) .* W_m_d(idx_medium).', 2);
    F_large  = sum(Q_export(:, idx_large)  .* W_m_d(idx_large).', 2);
else
    F_small  = sum(Q_export(:, idx_small)  .* V_m3(idx_small).'  .* W_m_d(idx_small).', 2);
    F_medium = sum(Q_export(:, idx_medium) .* V_m3(idx_medium).' .* W_m_d(idx_medium).', 2);
    F_large  = sum(Q_export(:, idx_large)  .* V_m3(idx_large).'  .* W_m_d(idx_large).', 2);
end

F_total_safe = max(F_total, 1e-30);
frac_small   = F_small  ./ F_total_safe;
frac_medium  = F_medium ./ F_total_safe;
frac_large   = F_large  ./ F_total_safe;

fig = figure('Name',[title_prefix ' Flux Analysis'], ...
    'Color','w','Position',[50 50 1400 900]);
turn_off_axes_toolbar(fig);

tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

% 1) epsilon
nexttile;
semilogy(t, eps_t, 'k-', 'LineWidth',1.8); hold on;
yline(1e-6, 'r--', 'LineWidth',1.2);
ylabel('\epsilon (W kg^{-1})');
title(title_prefix,'Interpreter','none');
grid on; box on; xlim([t(1) t(end)]);

% 2) relative export
nexttile;
plot(t, F_relative, 'b-', 'LineWidth',2.4); hold on;
yline(1,'k--','LineWidth',1.2);
ylabel('Relative export (F/F_{ref})');
title('Relative Export Flux','Interpreter','none');
grid on; box on; xlim([t(1) t(end)]);

% 3) fractions
nexttile;
plot(t, frac_large,  'LineWidth',2.2,'DisplayName','>2000 \mum'); hold on;
plot(t, frac_medium, 'LineWidth',2.2,'DisplayName','500–2000 \mum');
plot(t, frac_small,  'LineWidth',2.2,'DisplayName','<500 \mum');
ylim([0 1]);
ylabel('Export fraction');
xlabel('Time (days)');
legend('Location','eastoutside');
grid on; box on; xlim([t(1) t(end)]);

end

%% =================================================================
% EOF ANALYSIS (white at zero + larger layout + stable export)
%% =================================================================
function [fig_maps, fig_3d] = analyze_model_eofs_save( ...
    title_prefix, t, Y, D_um, cfg, t_forcing, eps_forcing, grid_obj)

t = t(:);
Nt = numel(t);
Ns = cfg.n_sections;
Nz = cfg.getNumLayers();
depths = cfg.getZ(); depths = depths(:);

eps_t = interp1(t_forcing(:), eps_forcing(:), t, 'previous', 'extrap');
eps_t = max(eps_t, 1e-14);

% reshape Nt x Nz x Ns
Y_3d = zeros(Nt, Nz, Ns);
for k = 1:Nz
    c1 = (k-1)*Ns + 1;
    c2 = k*Ns;
    Y_3d(:,k,:) = Y(:, c1:c2);
end

% EOF on log10(state)
Y_flat = reshape(Y_3d, Nt, Nz*Ns);
Y_log  = log10(max(Y_flat, 1e-20));
Y_mean = mean(Y_log, 1);
Y_anom = Y_log - Y_mean;

% remove dead columns
col_std = std(Y_anom, 0, 1);
keep = isfinite(col_std) & (col_std > 1e-12);
if ~any(keep), keep = true(1, size(Y_anom,2)); end
Y_anom2 = Y_anom(:, keep);

% PCA/SVD
try
    [coeffs2, scores2, latent] = pca(Y_anom2, 'Algorithm','svd','Economy',true,'Centered',false);
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
coeffs(keep,:) = coeffs2;

D_mm = D_um(:) / 1000;
logD = log10(D_mm);

mode1 = reshape(coeffs(:,1), Ns, Nz).';
mode2 = reshape(coeffs(:,2), Ns, Nz).';

% sign convention (optional)
if mean(mode1(1,end-3:end)) < 0
    mode1 = -mode1; scores(:,1) = -scores(:,1);
end
if mean(mode2(1,end-3:end)) < 0
    mode2 = -mode2; scores(:,2) = -scores(:,2);
end

% EOF colormap (white at 0) + symmetric scaling per mode
cmap = redblue_white0(256);
cl1  = max(abs(mode1(:))); if cl1 < 1e-12, cl1 = 1; end
cl2  = max(abs(mode2(:))); if cl2 < 1e-12, cl2 = 1; end

% PC time series
pc1 = scores(:,1); pc2 = scores(:,2);
pc1_n = pc1 ./ (max(abs(pc1)) + 1e-30);
pc2_n = pc2 ./ (max(abs(pc2)) + 1e-30);

% Export curve (normalized)
export_depth_idx = Nz;
c1 = (export_depth_idx-1)*Ns + 1;
c2 = export_depth_idx*Ns;
Q_export = Y(:, c1:c2);

r_cm    = grid_obj.getFractalRadii();
r_v     = grid_obj.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_obj.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;
V_m3    = (4/3) * pi * r_v.^3 * 1e-6;

if isprop(cfg,'state_is_biovolume') && logical(cfg.state_is_biovolume)
    F_total = sum(Q_export .* W_m_d.', 2);
else
    F_total = sum(Q_export .* V_m3.' .* W_m_d.', 2);
end
F_rel = F_total / max(max(F_total), 1e-30);

fig_maps = figure('Name',[title_prefix ' EOF Analysis'], ...
    'Color','w','Position',[50 50 1500 900]);
turn_off_axes_toolbar(fig_maps);
tiledlayout(2,2,'Padding','compact','TileSpacing','compact');

% forcing
nexttile;
semilogy(t, eps_t, 'k-','LineWidth',1.6);
ylabel('\epsilon (W kg^{-1})'); xlabel('Time (days)');
title('Turbulence forcing','Interpreter','none');
grid on; box on;

% PCs
nexttile;
plot(t, pc1_n, 'k-','LineWidth',1.5,'DisplayName','PC1'); hold on;
plot(t, pc2_n, 'r-','LineWidth',1.5,'DisplayName','PC2');
plot(t, F_rel,  'k--','LineWidth',1.2,'DisplayName','Export (norm)');
yline(0,'k:');
ylabel('Normalized'); xlabel('Time (days)');
title(sprintf('PC1 & PC2 (%.1f%% / %.1f%% var)', explained(1), explained(2)));
legend('Location','best'); grid on; box on;

% Mode 1
nexttile;
imagesc(logD, depths, mode1);
set(gca,'YDir','reverse');
colormap(gca, cmap);
caxis([-cl1 cl1]);
colorbar;
title(sprintf('Mode 1 (%.1f%% var)', explained(1)));
xlabel('D (mm), log_{10}'); ylabel('Depth (m)');
set(gca,'XTick',log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});
grid on;

% Mode 2
nexttile;
imagesc(logD, depths, mode2);
set(gca,'YDir','reverse');
colormap(gca, cmap);
caxis([-cl2 cl2]);
colorbar;
title(sprintf('Mode 2 (%.1f%% var)', explained(2)));
xlabel('D (mm), log_{10}'); ylabel('Depth (m)');
set(gca,'XTick',log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});
grid on;

% 3D slices: Full + Mode1 + Mode2
PSD_full_3d = reshape(Y_log, Nt, Nz, Ns);

Y1_anom = scores(:,1) * coeffs(:,1).';
Y1_log  = Y1_anom + Y_mean;
Y1_3d   = reshape(Y1_log, Nt, Nz, Ns);

Y2_anom = scores(:,2) * coeffs(:,2).';
Y2_log  = Y2_anom + Y_mean;
Y2_3d   = reshape(Y2_log, Nt, Nz, Ns);

clim_global = robust_clim(Y_log);

fig_3d = figure('Name',[title_prefix ' 3D slices'], ...
    'Color','w','Position',[50 50 1600 900]);
turn_off_axes_toolbar(fig_3d);
tiledlayout(3,1,'Padding','compact','TileSpacing','compact');

ax1 = nexttile; plot_psd_slice_subplot(ax1, t, depths, D_mm, PSD_full_3d, 'Full model (log_{10})', clim_global);
ax2 = nexttile; plot_psd_slice_subplot(ax2, t, depths, D_mm, Y1_3d,      'Mode 1 reconstruction', clim_global);
ax3 = nexttile; plot_psd_slice_subplot(ax3, t, depths, D_mm, Y2_3d,      'Mode 2 reconstruction', clim_global);

end

%% =================================================================
% Size spectra plots (log-log) at selected depths/times
% NOTE: For biovolume-state, it plots the state directly (shape is the point).
%% =================================================================
function figS = plot_size_spectra_selected(result, cfg, title_str, depths_m, times_d)

t = result.time(:);
Y = result.concentrations;

Ns = cfg.n_sections;
Nz = cfg.getNumLayers();
z  = cfg.getZ(); z = z(:);

% IMPORTANT: do NOT name this variable "grid"
gridObj = cfg.derive();
r_cm    = gridObj.getFractalRadii();
D_um    = 2 * r_cm * 1e4;

% depth indices
k_list = zeros(size(depths_m));
for i = 1:numel(depths_m)
    [~,k_list(i)] = min(abs(z - depths_m(i)));
end

% time indices
it_list = zeros(size(times_d));
for i = 1:numel(times_d)
    [~,it_list(i)] = min(abs(t - times_d(i)));
end

figS = figure('Color','w','Position',[50 50 1600 500],'Name',title_str);
turn_off_axes_toolbar(figS);
tiledlayout(1,numel(k_list),'Padding','compact','TileSpacing','compact');

for ii = 1:numel(k_list)
    k = k_list(ii);
    c1 = (k-1)*Ns + 1;
    c2 = k*Ns;

    nexttile; hold on;

    for jj = 1:numel(it_list)
        it = it_list(jj);
        q  = Y(it, c1:c2).';
        loglog(D_um, max(q,1e-30), 'LineWidth',1.6);
    end

    grid on; box on;
    xlabel('D [\mum]'); ylabel('State');
    title(sprintf('z = %.1f m', z(k)));
end

labels = arrayfun(@(x) sprintf('t=%.1f d', x), t(it_list), 'UniformOutput', false);
legend(labels, 'Location','eastoutside');
sgtitle(title_str,'Interpreter','none');

end

%% =================================================================
% 3D slice plotting helper
%% =================================================================
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
cb.Label.String = 'log_{10}';

xlabel('Time (days)');
ylabel('D (mm), log_{10}');
zlabel('Depth (m)');
title(title_str,'Interpreter','none');

set(gca,'YTick',log10([0.1 1 10]), 'YTickLabel',{'0.1','1','10'});
pbaspect([35 1 3]);
view([-60 6]);
axis tight; box on; grid on;
end

%% =================================================================
% Robust color limits helper
%% =================================================================
function clim = robust_clim(vals_in)
vals = vals_in(:);
vals = vals(isfinite(vals));
if isempty(vals)
    clim = [-20 -5];
    return;
end
vals = sort(vals);
n  = numel(vals);
lo = vals(max(1, round(0.02*n)));
hi = vals(round(0.98*n));
if hi <= lo, hi = lo + 1; end
clim = [lo hi];
end

%% =================================================================
% Load epsilon forcing helper
%% =================================================================
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

%% =================================================================
% Diverging colormap: blue -> white -> red (white at 0)
%% =================================================================
function cmap = redblue_white0(n)
if nargin<1, n = 256; end
n2 = floor(n/2);

% blue -> white
r1 = linspace(0,1,n2)';  
g1 = linspace(0,1,n2)';  
b1 = ones(n2,1);

% white -> red
r2 = ones(n-n2,1);
g2 = linspace(1,0,n-n2)'; 
b2 = linspace(1,0,n-n2)';

cmap = [r1 g1 b1; r2 g2 b2];
end

%% =================================================================
% Turn off toolbar so exportgraphics/save does not warn
%% =================================================================
function turn_off_axes_toolbar(fig)
% remove figure toolbar + axes toolbars (exportgraphics warning fix)

try, fig.ToolBar = 'none'; end
try, fig.MenuBar = 'none'; end

ax = findall(fig,'type','axes');
for k = 1:numel(ax)
    try
        ax(k).Toolbar.Visible = 'off';
    catch
    end
end
end