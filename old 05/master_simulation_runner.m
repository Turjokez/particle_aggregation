function master_simulation_runner
% MASTER_SIMULATION_RUNNER  – Full tuned run
%
%  Runs:
%    B1) Atlantic column with observed epsilon
%    B2) Square-wave turbulence test
%    B3) EOF analysis (attenuation + seesaw) from Atlantic run
%
%  Tuned parameters (based on Adrian feedback):
%    - attenuation_rate (mu)      = 0.12 day^-1
%    - attenuation_depth_factor   = 0.5   (mu increases with depth)
%    - NPP_rate                   = 3e-5  day^-1
%    - NPP finite bloom: on until day 10, off afterwards
%    - Stronger fragmentation / disaggregation

clc;
mu_tuned        = 0.12;
NPP_tuned       = 3e-5;
NPP_stop_day    = 10;
depth_slope     = 0.5;

outDir = 'fig_modelMain_tuned';
if ~exist(outDir,'dir')
    mkdir(outDir);
end

fprintf('======================================================\n');
fprintf('   PARTICLE AGGREGATION MODEL - FULL RUN (TUNED)\n');
fprintf('   mu = %.2f day^-1, NPP = %.1e (finite bloom)\n', mu_tuned, NPP_tuned);
fprintf('======================================================\n\n');

%% ---------------------------------------------------------------
%  B1: Atlantic column (observed epsilon)
% ---------------------------------------------------------------
fprintf('--- B1: Atlantic column (tuned mu/NPP) ---\n');

[simA, cfgA, eps_runA] = run_atlantic_column(mu_tuned, NPP_tuned, ...
    NPP_stop_day, depth_slope);

figA = make_flux_plots(simA.result, simA.grid, cfgA, eps_runA, ...
    sprintf('Atlantic Observed (mu=%.2f, NPP=%.1e)', mu_tuned, NPP_tuned));

saveas(figA, fullfile(outDir, 'B1_Atlantic_flux.png'));

%% ---------------------------------------------------------------
%  B2: Square-wave turbulence test
% ---------------------------------------------------------------
fprintf('\n--- B2: Square-wave turbulence test (tuned mu/NPP) ---\n');

[simB, cfgB, eps_runB] = run_square_wave(mu_tuned, NPP_tuned, ...
    NPP_stop_day, depth_slope);

figB = make_flux_plots(simB.result, simB.grid, cfgB, eps_runB, ...
    sprintf('Square Wave Test (mu=%.2f, NPP=%.1e)', mu_tuned, NPP_tuned));

saveas(figB, fullfile(outDir, 'B2_squareWave_flux.png'));

%% ---------------------------------------------------------------
%  B3: EOF analysis using Atlantic run
% ---------------------------------------------------------------
fprintf('\n--- B3: EOF analysis (Atlantic, tuned mu/NPP) ---\n');

tA   = simA.result.time;
YA   = simA.result.concentrations;
D_um = 2 * simA.grid.getFractalRadii() * 1e4;

[figC_maps, figC_3d] = analyze_model_eofs_save( ...
    'Atlantic', tA, YA, D_um, cfgA, eps_runA, outDir, 'B3_EOF');

saveas(figC_maps, fullfile(outDir, 'B3_EOF_maps_PCs.png'));
saveas(figC_3d,   fullfile(outDir, 'B3_EOF_3D_slices.png'));

fprintf('\nAll main figures saved in folder: %s\n', outDir);
fprintf('Done.\n');
end

%% ========================================================================
%  B1: Atlantic column run
%  ========================================================================
function [sim, cfg, eps_run] = run_atlantic_column(mu_val, NPP_val, ...
    NPP_stop_day, depth_slope)

if ~isfile('epsilon_daily.mat')
    error('Missing epsilon_daily.mat in current folder.');
end

S = load('epsilon_daily.mat');
[t_days, eps_raw] = grab_data_simple(S);

% first 60 days (eddy-core window)
mask    = t_days <= 60;
t_run   = t_days(mask);
eps_run = eps_raw(mask);

fprintf('Running Atlantic simulation for %.1f days (N=%d steps)...\n', ...
    range(t_run), numel(t_run));

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 200;
cfg.dz         = 10;

% stronger nonlinear disaggregation / fragmentation
cfg.disagg_use_nonlinear  = true;
cfg.disagg_kmax_a         = 0.99;
cfg.disagg_beta           = 1.0;
cfg.disagg_redistribute_p = 2.0;
cfg.c3                    = 0.2;   % stronger linear breakup
cfg.c4                    = 1.6;

% tuned attenuation & NPP (finite bloom)
cfg.attenuation_rate         = mu_val;
cfg.attenuation_depth_factor = depth_slope;

cfg.use_NPP        = true;
cfg.NPP_rate       = NPP_val;
cfg.NPP_t_step     = NPP_stop_day;  % after this time...
cfg.NPP_rate_after = 0;             % ...turn NPP off

sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_run, eps_run);

sim.run('tspan', t_run);
end

%% ========================================================================
%  B2: Square-wave turbulence run
%  ========================================================================
function [sim, cfg, eps_run] = run_square_wave(mu_val, NPP_val, ...
    NPP_stop_day, depth_slope)

T_total = 15;      % days
dt      = 0.05;
t_run   = (0:dt:T_total).';

% low background epsilon
eps_run = 1e-9 * ones(size(t_run));

% three pulses (~1 day each)
pulse_edges = [0.5  1.5;
               5.0  6.0;
               9.5 10.5];
for k = 1:size(pulse_edges,1)
    mask = (t_run >= pulse_edges(k,1)) & (t_run < pulse_edges(k,2));
    eps_run(mask) = 5e-6;
end

fprintf('Running Square-wave simulation for %.1f days (N=%d steps)...\n', ...
    range(t_run), numel(t_run));

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 200;
cfg.dz         = 10;

% same tuned breakup settings
cfg.disagg_use_nonlinear  = true;
cfg.disagg_kmax_a         = 0.99;
cfg.disagg_beta           = 1.0;
cfg.disagg_redistribute_p = 2.0;
cfg.c3                    = 0.2;
cfg.c4                    = 1.6;

cfg.attenuation_rate         = mu_val;
cfg.attenuation_depth_factor = depth_slope;

cfg.use_NPP        = true;
cfg.NPP_rate       = NPP_val;
cfg.NPP_t_step     = NPP_stop_day;
cfg.NPP_rate_after = 0;

sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_run, eps_run);
sim.run('tspan', t_run);
end

%% ========================================================================
%  Shared FLUX + SIZE-CLASS plotting
%  ========================================================================
function fig = make_flux_plots(result, grid_obj, cfg, eps_run, title_prefix)

t  = result.time(:);
Y  = result.concentrations;
Nt = numel(t);
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;

% grid / settling
r_cm    = grid_obj.getFractalRadii();
r_v     = grid_obj.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_obj.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;              % cm/s -> m/day
V_m3    = (4/3) * pi * r_v.^3 * 1e-6;           % cm^3 -> m^3
D_um    = 2 * r_cm * 1e4;                       % µm

% export flux at bottom layer
export_depth_idx = Nz;
Q_export = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);
F_total  = sum(Q_export .* V_m3.' .* W_m_d.', 2);   % Nt x 1

% steady-state flux from low-epsilon periods
F_ss_idx = find(eps_run < 1e-10, 10);
if isempty(F_ss_idx)
    F_ss_idx = 1:min(10,Nt);
end
F_ss = mean(F_total(F_ss_idx));
if F_ss == 0
    F_ss = 1;
end
F_relative = F_total / F_ss;

% size classes
D_small_max  = 200;    % <200 µm
D_medium_max = 1000;   % 200–1000 µm
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

% panel 1: epsilon
subplot(3,1,1);
semilogy(t, max(eps_run,1e-12), 'k-', 'LineWidth',1.2); hold on;
plot(t([1 end]), [1e-6 1e-6], 'r--', 'LineWidth',1);
ylabel('\epsilon (W kg^{-1})');
title(title_prefix);
axis tight; grid on;

% panel 2: relative flux
subplot(3,1,2);
plot(t, F_relative, 'b-', 'LineWidth',2); hold on;
yline(1,'k--','LineWidth',1);
ylabel('Relative Flux (F/F_{ss})');
title('Relative Export Flux');
axis tight; grid on;

% panel 3: flux fractions
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

%% ========================================================================
%  EOF analysis + 3D slices (saves two figures)
%  ========================================================================
function [fig_maps, fig_3d] = analyze_model_eofs_save( ...
    title_prefix, t, Y, D_um, cfg, eps_run, outDir, tag)

Nt = numel(t);
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;
depths = (1:Nz) * cfg.dz;

% reshape [time, depth, size]
Y_3d = zeros(Nt, Nz, Ns);
for k = 1:Nz
    idx_start = (k-1)*Ns + 1;
    idx_end   = k*Ns;
    Y_3d(:,k,:) = Y(:, idx_start:idx_end);
end

% log10 PSD
Y_flat = reshape(Y_3d, Nt, Nz*Ns);
Y_log  = log10(max(Y_flat, 1e-20));
Y_mean = mean(Y_log, 1);
Y_anom = Y_log - Y_mean;

% PCA
[coeffs, scores, ~, ~, explained] = pca(Y_anom);

D_mm = D_um(:) / 1000;
logD = log10(D_mm);

mode1 = reshape(coeffs(:,1), Ns, Nz).';
mode2 = reshape(coeffs(:,2), Ns, Nz).';

% sign conventions
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

% ---- Figure: maps + PCs ----
fig_maps = figure('Name',[title_prefix ' EOF Analysis'], ...
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

% normalized flux
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;
export_depth_idx = Nz;
Q_export = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);

sim_tmp  = CoagulationSimulation(cfg);
r_cm    = sim_tmp.grid.getFractalRadii();
r_v     = sim_tmp.grid.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, sim_tmp.grid.setcon);
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
title(sprintf('Mode 1 (Attenuation, %.1f%%%% var)', explained(1)));
xlabel('D (mm), log_{10} scale');
ylabel('Depth (m)');
set(gca,'XTick',log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});

subplot(2,2,4);
imagesc(logD, depths, mode2);
set(gca,'YDir','reverse');
colormap(gca, turbo);
cb = colorbar; cb.Label.String = 'EOF2 loading';
title(sprintf('Mode 2 (Seesaw, %.1f%%%% var)', explained(2)));
xlabel('D (mm), log_{10} scale');
ylabel('Depth (m)');
set(gca,'XTick',log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});

% ---- Figure: 3D slices ----
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

fig_3d = figure('Name',[title_prefix ' 3D Reconstruction'], ...
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
%  Helper: 3D slice "fence" plot
%  ========================================================================
function plot_psd_slice_subplot(ax_handle, t_days, depths, D_mm, PSD_3d, ...
    title_str, c_limits)

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
%  Utility: load epsilon file robustly
%  ========================================================================
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