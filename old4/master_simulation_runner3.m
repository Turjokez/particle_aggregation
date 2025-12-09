function master_simulation_runner3
% MASTER_SIMULATION_RUNNER3
%
% Version that uses tuned attenuation (MU_STAR) and NPP (NPP_STAR)
% in all three experiments:
%   [1] Real Atlantic Data (Flux + Size Classes)
%   [2] Square Wave Test (Idealized Turbulence Pulses)
%   [3] EOF Analysis (Attenuation + Seesaw Modes)
%
% Turjo – change MU_STAR and NPP_STAR here when you decide on the
% final tuned pair (ideally consistent with CASE 3 in runner2).

% -------------------------------------------------------------------------
% Tuned parameters (edit these two numbers if you want to retune)
% -------------------------------------------------------------------------
MU_STAR  = 0.10;    % attenuation rate [d^-1]
NPP_STAR = 5e-5;    % NPP source rate (upper layer, smallest bin)

clc;
fprintf('======================================================\n');
fprintf('   PARTICLE AGGREGATION MODEL - RESEARCH (TUNED)     \n');
fprintf('======================================================\n');
fprintf('Using: mu = %.3f d^{-1}, NPP = %.2e\n', MU_STAR, NPP_STAR);
fprintf('======================================================\n');
fprintf('Select an experiment:\n');
fprintf('  [1] Real Atlantic Data (Flux + Size Classes)\n');
fprintf('  [2] Square Wave Test (Idealized Turbulence Pulses)\n');
fprintf('  [3] EOF Analysis (Attenuation + Seesaw Modes)\n');
fprintf('======================================================\n');

choice = input('Enter number [1-3]: ');

switch choice
    case 1
        run_column_experiment_tuned(MU_STAR, NPP_STAR);
    case 2
        run_square_wave_experiment_tuned(MU_STAR, NPP_STAR);
    case 3
        run_eof_analysis_experiment_tuned(MU_STAR, NPP_STAR);
    otherwise
        fprintf('Invalid selection.\n');
end
end

%% ========================================================================
%  EXPERIMENT 1: REAL ATLANTIC DATA – FLUX + SIZE CLASSES
% ========================================================================
function run_column_experiment_tuned(MU_STAR, NPP_STAR)
fprintf('\n=== 1-D COLUMN SIMULATION [ATLANTIC, TUNED MU/NPP] ===\n');

if ~isfile('epsilon_daily.mat')
    error('Please ensure epsilon_daily.mat is in the folder.');
end
S = load('epsilon_daily.mat');
[t_days, eps_raw] = grab_data_simple(S);

% Eddy-core window (same as before; you can adjust if needed)
mask    = t_days <= 27;
t_run   = t_days(mask);
eps_run = eps_raw(mask);

fprintf('Running simulation for %.1f days (N=%d steps)...\n', ...
    range(t_run), numel(t_run));

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 200;
cfg.dz         = 10;

cfg.disagg_use_nonlinear  = true;
cfg.disagg_kmax_a         = 0.95;
cfg.disagg_beta           = 1.0;
cfg.disagg_redistribute_p = 1.5;

cfg.attenuation_rate = MU_STAR;
cfg.use_NPP          = true;
cfg.NPP_rate         = NPP_STAR;

cfg.t_init  = t_run(1);
cfg.t_final = t_run(end);
cfg.delta_t = t_run(2) - t_run(1);

sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_run, eps_run);

sim.run('tspan', t_run);

make_flux_plots(sim.result, sim.grid, cfg, eps_run, ...
    sprintf('Atlantic Observed (mu=%.3f, NPP=%.1e)', MU_STAR, NPP_STAR));
end

%% ========================================================================
%  EXPERIMENT 2: IDEALIZED SQUARE-WAVE TURBULENCE + FLUX
% ========================================================================
function run_square_wave_experiment_tuned(MU_STAR, NPP_STAR)
fprintf('\n=== SQUARE WAVE EXPERIMENT (TUNED MU/NPP) ===\n');

T_total = 15;               % days
dt      = 0.05;
t_run   = (0:dt:T_total).';

% Low background epsilon
eps_run = 1e-9 * ones(size(t_run));

% Three ~1-day pulses
pulse_edges = [0.5 1.5;
               5.0 6.0;
               9.5 10.5];
for k = 1:size(pulse_edges,1)
    mask = (t_run >= pulse_edges(k,1)) & (t_run < pulse_edges(k,2));
    eps_run(mask) = 5e-6;
end

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 200;
cfg.dz         = 10;

cfg.disagg_use_nonlinear  = true;
cfg.disagg_kmax_a         = 0.95;
cfg.disagg_beta           = 1.0;
cfg.disagg_redistribute_p = 1.5;

cfg.attenuation_rate = MU_STAR;
cfg.use_NPP          = true;
cfg.NPP_rate         = NPP_STAR;

cfg.t_init  = t_run(1);
cfg.t_final = t_run(end);
cfg.delta_t = t_run(2) - t_run(1);

sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_run, eps_run);
sim.run('tspan', t_run);

make_flux_plots(sim.result, sim.grid, cfg, eps_run, ...
    sprintf('Square Wave Test (mu=%.3f, NPP=%.1e)', MU_STAR, NPP_STAR));
end

%% ========================================================================
%  EXPERIMENT 3: EOF ANALYSIS + 3D PSD SLICES (TUNED)
% ========================================================================
function run_eof_analysis_experiment_tuned(MU_STAR, NPP_STAR)
fprintf('\n=== MODEL EOF ANALYSIS [TUNED MU/NPP] ===\n');

[out, cfg, eps_run] = run_simulation_internal_tuned(MU_STAR, NPP_STAR);

t    = out.t_days;
Y    = out.concentrations;
D_um = out.D_um;

analyze_model_eofs('Atlantic (tuned)', t, Y, D_um, cfg, eps_run);
end

%% ========================================================================
%  INTERNAL RUN FOR EOF ANALYSIS (Atlantic forcing + tuned mu/NPP)
% ========================================================================
function [out, cfg, eps_run] = run_simulation_internal_tuned(MU_STAR, NPP_STAR)

if ~isfile('epsilon_daily.mat')
    error('Missing epsilon_daily.mat for internal EOF run.');
end
S = load('epsilon_daily.mat');
[t_days, eps_raw] = grab_data_simple(S);

mask    = t_days <= 27;
t_run   = t_days(mask);
eps_run = eps_raw(mask);

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 200;
cfg.dz         = 10;

cfg.disagg_use_nonlinear  = true;
cfg.disagg_kmax_a         = 0.95;
cfg.disagg_beta           = 1.0;
cfg.disagg_redistribute_p = 1.5;

cfg.attenuation_rate = MU_STAR;
cfg.use_NPP          = true;
cfg.NPP_rate         = NPP_STAR;

cfg.t_init  = t_run(1);
cfg.t_final = t_run(end);
cfg.delta_t = t_run(2) - t_run(1);

sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_run, eps_run);
sim.run('tspan', t_run);

out_temp           = sim.exportMinimalOutputs();
out.t_days         = sim.result.time;
out.concentrations = sim.result.concentrations;
out.D_um           = out_temp.D_um;
end

%% ========================================================================
%  FLUX + SIZE-CLASS PLOTTING (shared by exp 1 & 2)
% ========================================================================
function make_flux_plots(result, grid_obj, cfg, eps_run, title_prefix)

t  = result.time(:);
Y  = result.concentrations;
Nt = numel(t);
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;

r_cm    = grid_obj.getFractalRadii();
r_v     = grid_obj.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_obj.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;
V_m3    = (4/3) * pi * r_v.^3 * 1e-6;
D_um    = 2 * r_cm * 1e4;

export_depth_idx = Nz;
Q_export = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);

F_total = sum(Q_export .* V_m3.' .* W_m_d.', 2);

F_ss_idx = find(eps_run < 1e-10, 10);
if isempty(F_ss_idx)
    F_ss_idx = 1:min(10, Nt);
end
F_ss = mean(F_total(F_ss_idx));
if F_ss == 0, F_ss = 1; end
F_relative = F_total / F_ss;

D_small_max  = 200;
D_medium_max = 1000;
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

figure('Name',[title_prefix ' – Flux Analysis'], ...
       'Color','w','Position',[100 100 800 900]);

subplot(3,1,1);
semilogy(t, max(eps_run,1e-12), 'k-', 'LineWidth',1.2); hold on;
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
plot(t, frac_large,  'r-','LineWidth',2.2,'DisplayName','Large');
hold on;
plot(t, frac_medium,'g-','LineWidth',2.2,'DisplayName','Medium');
plot(t, frac_small, 'b-','LineWidth',2.2,'DisplayName','Small');
ylim([0 1]);
ylabel('Flux Fraction');
xlabel('Time (days)');
legend('Location','northeast');
grid on; box on;
end

%% ========================================================================
%  EOF ANALYSIS + 3D SLICE PLOTS
% ========================================================================
function analyze_model_eofs(title_prefix, t, Y, D_um, cfg, eps_run)

Nt = numel(t);
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;
depths = (1:Nz) * cfg.dz;

Y_3d = zeros(Nt, Nz, Ns);
for k = 1:Nz
    idx_start = (k-1)*Ns + 1;
    idx_end   = k*Ns;
    Y_3d(:,k,:) = Y(:, idx_start:idx_end);
end

Y_flat = reshape(Y_3d, Nt, Nz*Ns);
Y_log  = log10(max(Y_flat, 1e-20));
Y_mean = mean(Y_log, 1);
Y_anom = Y_log - Y_mean;

[coeffs, scores, ~, ~, explained] = pca(Y_anom);

D_mm = D_um(:) / 1000;
logD = log10(D_mm);

mode1 = reshape(coeffs(:,1), Ns, Nz).';
mode2 = reshape(coeffs(:,2), Ns, Nz).';

if mean(mode1(1, end-3:end)) < 0
    mode1       = -mode1;
    coeffs(:,1) = -coeffs(:,1);
    scores(:,1) = -scores(:,1);
end
if mean(mode2(1, end-3:end)) < 0
    mode2       = -mode2;
    coeffs(:,2) = -coeffs(:,2);
    scores(:,2) = -scores(:,2);
end

figure('Name',[title_prefix ' EOF Analysis (tuned)'], ...
       'Color','w','Position',[50 50 1200 850]);

subplot(2,2,1);
semilogy(t, max(eps_run,1e-13), 'k-','LineWidth',1.2);
ylabel('\epsilon (W kg^{-1})');
xlabel('Time (days)');
title('Turbulence Forcing');
axis tight; grid on;

subplot(2,2,2);
pc1 = scores(:,1);
pc2 = scores(:,2);
pc1_n = pc1 ./ max(abs(pc1));
pc2_n = pc2 ./ max(abs(pc2));

% flux (bottom layer) for reference
export_depth_idx = Nz;
Q_export = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);

grid_tmp = cfg.derive();
r_cm  = grid_tmp.getFractalRadii();
r_v   = grid_tmp.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_tmp.setcon);
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
title(sprintf('PC1 & PC2 (%.1f%% / %.1f%% var)', ...
    explained(1), explained(2)));
legend('Location','best');
axis tight; grid on;

subplot(2,2,3);
imagesc(logD, depths, mode1);
set(gca,'YDir','reverse');
colormap(gca, turbo);
cb = colorbar;
cb.Label.String = 'EOF1 loading';
title(sprintf('Mode 1 (Attenuation, %.1f%% var)', explained(1)));
xlabel('D (mm), log_{10} scale');
ylabel('Depth (m)');
set(gca,'XTick',log10([0.1 1 10]), ...
        'XTickLabel',{'0.1','1','10'});

subplot(2,2,4);
imagesc(logD, depths, mode2);
set(gca,'YDir','reverse');
colormap(gca, turbo);
cb = colorbar;
cb.Label.String = 'EOF2 loading';
title(sprintf('Mode 2 (Seesaw, %.1f%% var)', explained(2)));
xlabel('D (mm), log_{10} scale');
ylabel('Depth (m)');
set(gca,'XTick',log10([0.1 1 10]), ...
        'XTickLabel',{'0.1','1','10'});

% ---- Reconstruction plots (full, mode1, mode2) --------------------------
PSD_full_log = Y_log;
PSD_full_3d  = reshape(PSD_full_log, Nt, Nz, Ns);

Y1_anom = scores(:,1) * coeffs(:,1).';
Y1_log  = Y1_anom + Y_mean;
Y1_3d   = reshape(Y1_log, Nt, Nz, Ns);

Y2_anom = scores(:,2) * coeffs(:,2).';
Y2_log  = Y2_anom + Y_mean;
Y2_3d   = reshape(Y2_log, Nt, Nz, Ns);

vals = sort(PSD_full_log(:));
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

figure('Name',[title_prefix ' 3D Reconstruction (tuned)'], ...
       'Color','w','Position',[50 50 1100 950]);

ax1 = subplot(3,1,1);
plot_psd_slice_subplot(ax1, t, depths, D_mm, PSD_full_3d, ...
    'Full  Model (log_{10} Part Vol)', clim_global);

ax2 = subplot(3,1,2);
plot_psd_slice_subplot(ax2, t, depths, D_mm, Y1_3d, ...
    'Mode 1 reconstruction (attenuation)', clim_global);

ax3 = subplot(3,1,3);
plot_psd_slice_subplot(ax3, t, depths, D_mm, Y2_3d, ...
    'Mode 2 reconstruction (seesaw)', clim_global);
end

%% ========================================================================
function plot_psd_slice_subplot(ax_handle, t_days, depths, D_mm, PSD_3d, ...
    title_str, c_limits)

axes(ax_handle); %#ok<LAXES>
hold on;

t_days = t_days(:);
depths = depths(:);
logD   = log10(D_mm(:));

V = permute(PSD_3d, [3 1 2]);
[X,Y,Z] = meshgrid(t_days, logD, depths);

t_min = ceil(min(t_days));
t_max = floor(max(t_days));
xs    = t_min:1:t_max;
xslice = fliplr(xs);
yslice = [];
zslice = [];

h = slice(X,Y,Z,V,xslice,yslice,zslice);
set(h,'EdgeColor','none','FaceAlpha',0.95);

set(gca,'ZDir','reverse');
colormap(gca, turbo);

if numel(c_limits)==2 && c_limits(2) > c_limits(1)
    caxis(c_limits);
end

cb = colorbar('Location','eastoutside');
cb.Label.String = 'log_{10} Part Vol (ppmV mm^{-1})';

xlabel('Time (days)');
ylabel('D (mm), log_{10} scale');
zlabel('Depth (m)');
title(title_str);

set(gca,'YTick',log10([0.1 1 10]), ...
        'YTickLabel',{'0.1','1','10'});

pbaspect([35 1 3]);
view([-60 5]);
axis tight; box on; grid on;
end

%% ========================================================================
%  SMALL UTILITY: LOAD EPSILON FILES ROBUSTLY
% ========================================================================
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

t = t(:) - t(1);
n = min(numel(t), numel(e));
t = t(1:n);
e = e(1:n);
e = max(e, 1e-12);
end