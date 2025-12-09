function master_simulation_runner
% MASTER_SIMULATION_RUNNER
% All-in-one script for Particle Aggregation Project.
%
% Version: Atlantic forcing + square-wave test + EOF analysis
%          with attenuation and size-dependent disaggregation.
%
% Menu:
%   [1] Atlantic forcing  (flux + size classes)
%   [2] Square-wave test  (flux + size classes)
%   [3] EOF analysis on Atlantic run
%
% Note: Attenuation (mu) is tested over several values internally,
%       but EOF analysis uses mu = 0.05 day^-1 as the reference case.

clc;
fprintf('======================================================\n');
fprintf('   PARTICLE AGGREGATION MODEL - RESEARCH EDITION\n');
fprintf('======================================================\n');
fprintf('Select an experiment:\n');
fprintf('  [1] Real Atlantic Data (Flux + Size Classes)\n');
fprintf('  [2] Square Wave Test (Idealized Turbulence Forcing)\n');
fprintf('  [3] EOF Analysis (Atlantic Column with Attenuation)\n');
fprintf('======================================================\n');

choice = input('Enter number [1-3]: ');

switch choice
    case 1
        run_column_experiment();
    case 2
        run_square_wave_experiment();
    case 3
        analyze_model_eofs('Atlantic');
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

% Limit to first 60 days (or whatever you like)
mask    = t_days <= 60;
t_run   = t_days(mask);
eps_run = eps_raw(mask);

fprintf('Running simulation for %.1f days (N = %d steps)...\n', ...
    range(t_run), numel(t_run));

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 200;
cfg.dz         = 10;

% --- size-dependent disaggregation (nonlinear) ---------------------------
cfg.disagg_use_nonlinear = true;
cfg.disagg_kmax_a        = 0.95;
cfg.disagg_beta          = 1.0;
cfg.disagg_redistribute_p = 1.5;

% --- attenuation term (mu) ----------------------------------------------
cfg.attenuation_rate = 0.05;       % day^-1 (reference)

% --- simple constant NPP source (as before) ------------------------------
cfg.use_NPP  = true;
cfg.NPP_rate = 2e-4;

sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_run, eps_run);

tic;
sim.run('tspan', t_run);
fprintf('Simulation finished in %.2f seconds.\n', toc);

make_comprehensive_plots(sim.result, sim.grid, cfg, eps_run, ...
    'Atlantic Observed');
end

% =========================================================================
% === EXPERIMENT 2: SQUARE WAVE TEST ======================================
% =========================================================================
function run_square_wave_experiment
fprintf('\n=== SQUARE WAVE EXPERIMENT (WITH ATTENUATION) ===\n');

T_total = 15;        % days
dt      = 0.05;      % days
t_run   = (0:dt:T_total).';

% Background epsilon
eps_run = 1e-9 * ones(size(t_run));

% Three square pulses (~1 day wide)
pulse_edges = [0.5 1.5;
               5.0 6.0;
               9.5 10.5];

for k = 1:size(pulse_edges, 1)
    mask = (t_run >= pulse_edges(k,1)) & (t_run < pulse_edges(k,2));
    eps_run(mask) = 5e-6;
end

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 200;
cfg.dz         = 10;

% size-dependent disaggregation
cfg.disagg_use_nonlinear = true;
cfg.disagg_kmax_a        = 0.95;
cfg.disagg_beta          = 1.0;
cfg.disagg_redistribute_p = 1.5;

% attenuation
cfg.attenuation_rate = 0.05;   % same reference mu

cfg.use_NPP  = true;
cfg.NPP_rate = 2e-4;

sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_run, eps_run);
sim.run('tspan', t_run);

make_comprehensive_plots(sim.result, sim.grid, cfg, eps_run, ...
    'Square Wave Test');
end

% =========================================================================
% === CORE PLOTTING FUNCTION (FLUX & SIZE CLASSES) ========================
% =========================================================================
function make_comprehensive_plots(result, sim_grid_obj, cfg, eps_run, ...
    title_prefix)

t  = result.time(:);           % Nt
Y  = result.concentrations;    % [Nt x (Nz*Ns)]
Nt = numel(t);
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;

% --- particle radii, volumes, settling speeds ---------------------------
r_cm = sim_grid_obj.getFractalRadii();
r_v  = sim_grid_obj.getConservedRadii();

ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, sim_grid_obj.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;          % cm/s -> m/day
V_m3    = (4/3) * pi * r_v.^3 * 1e-6;       % cm^3 -> m^3
D_um    = 2 * r_cm * 1e4;                   % microns

% --- export flux at bottom level ----------------------------------------
export_depth_idx = Nz;
Q_export = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);

F_total = sum(Q_export .* (V_m3.' .* W_m_d.'), 2);   % Nt x 1

% steady-state flux (low epsilon periods)
eps_floor = 1e-10;
idx_ss = find(eps_run < eps_floor);
if numel(idx_ss) < 5
    idx_ss = 1:min(10, Nt);
end
F_ss = mean(F_total(idx_ss));
if F_ss == 0
    F_ss = 1;
end

F_relative = F_total / F_ss;

% --- size classes --------------------------------------------------------
D_small_max  = 200;   % µm
D_medium_max = 1000;  % µm

idx_small  = (D_um < D_small_max);
idx_medium = (D_um >= D_small_max) & (D_um < D_medium_max);
idx_large  = (D_um >= D_medium_max);

F_small  = sum(Q_export(:, idx_small)  .* (V_m3(idx_small).'  .* W_m_d(idx_small).'), 2);
F_medium = sum(Q_export(:, idx_medium) .* (V_m3(idx_medium).' .* W_m_d(idx_medium).'), 2);
F_large  = sum(Q_export(:, idx_large)  .* (V_m3(idx_large).'  .* W_m_d(idx_large).'), 2);

F_total_safe = max(F_total, 1e-30);
frac_small   = F_small  ./ F_total_safe;
frac_medium  = F_medium ./ F_total_safe;
frac_large   = F_large  ./ F_total_safe;

% --- figure --------------------------------------------------------------
figure('Name', [title_prefix ' Flux Analysis'], ...
       'Color', 'w', 'Position', [100 100 800 900]);

% Panel 1: turbulence forcing
subplot(3,1,1);
semilogy(t, max(eps_run,1e-12), 'k-', 'LineWidth', 1.2); hold on;
yline(1e-6, 'r--', 'LineWidth', 1);
ylabel('\epsilon (W kg^{-1})');
title([title_prefix ' Simulation']);
xlim([t(1) t(end)]);
grid on; box on;

% Panel 2: relative total export flux
subplot(3,1,2);
plot(t, F_relative, 'b-', 'LineWidth', 2);
hold on;
yline(1, 'k--');
ylabel('Relative Flux (F/F_{ss})');
title('Relative Export Flux');
xlim([t(1) t(end)]);
grid on; box on;

% Panel 3: flux fractions
subplot(3,1,3);
plot(t, frac_large,  'r-',  'LineWidth', 2.2, 'DisplayName','Large');
hold on;
plot(t, frac_medium, 'g--', 'LineWidth', 2.2, 'DisplayName','Medium');
plot(t, frac_small,  'b:',  'LineWidth', 2.2, 'DisplayName','Small');
ylim([0 1]);
ylabel('Flux Fraction');
xlabel('Time (days)');
legend('Location','southwest');
xlim([t(1) t(end)]);
grid on; box on;
end

% =========================================================================
% === EOF ANALYSIS + 3D RECONSTRUCTIONS ===================================
% =========================================================================
function analyze_model_eofs(title_prefix, t, Y, cfg, eps_run, F_relative)

% If called with only the title, run an internal Atlantic simulation
if nargin < 2
    fprintf('\n=== MODEL EOF ANALYSIS [ATLANTIC, WITH ATTENUATION] ===\n');
    [out, cfg, eps_run, flux_struct] = run_simulation_internal();
    t          = out.t_days;
    Y          = out.concentrations;
    F_relative = flux_struct.F_relative;
    title_prefix = 'Atlantic';
end

t = t(:);
Nt = numel(t);

Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;
depths = (1:Nz) * cfg.dz;

% --- reshape to [Nt x Nz x Ns] ------------------------------------------
Y_3d = zeros(Nt, Nz, Ns);
for k = 1:Nz
    idx_start = (k-1)*Ns + 1;
    idx_end   = k*Ns;
    Y_3d(:,k,:) = Y(:, idx_start:idx_end);
end

% --- log10 PSD & anomalies ----------------------------------------------
Y_flat = reshape(Y_3d, Nt, Nz*Ns);             % [Nt x (Nz*Ns)]
Y_log  = log10(max(Y_flat, 1e-20));
Y_mean = mean(Y_log, 1);
Y_anom = Y_log - Y_mean;

[coeffs, scores, ~, ~, explained] = pca(Y_anom);     % coeffs: [Nz*Ns x Nmodes]

% --- size axis in mm -----------------------------------------------------
sim_tmp = CoagulationSimulation(cfg);   % only to get size grid
D_um    = 2 * sim_tmp.grid.getFractalRadii() * 1e4;
D_mm    = D_um / 1000;
logD    = log10(D_mm);

% --- first two modes reshaped to [Nz x Ns] ------------------------------
mode1 = reshape(coeffs(:,1), Ns, Nz).';
mode2 = reshape(coeffs(:,2), Ns, Nz).';

% Flip signs so large near-surface loadings are positive
if mean(mode1(1,end-2:end)) < 0
    mode1       = -mode1;
    coeffs(:,1) = -coeffs(:,1);
    scores(:,1) = -scores(:,1);
end
if mean(mode2(1,end-2:end)) < 0
    mode2       = -mode2;
    coeffs(:,2) = -coeffs(:,2);
    scores(:,2) = -scores(:,2);
end

% --- FIGURE 1: epsilon, PCs, mode maps ----------------------------------
figure('Name',[title_prefix ' EOF Analysis'], 'Color','w', ...
       'Position',[50 50 1200 850]);

% (a) turbulence forcing
subplot(2,2,1);
semilogy(t, max(eps_run,1e-12), 'k-', 'LineWidth',1.2);
ylabel('\epsilon (W kg^{-1})');
xlabel('Time (days)');
title('Turbulence Forcing');
xlim([t(1) t(end)]);
grid on; box on;

% (b) PC1, PC2, and normalized relative flux
subplot(2,2,2);
pc1 = scores(:,1);
pc2 = scores(:,2);

pc1_n = pc1 ./ max(abs(pc1));
pc2_n = pc2 ./ max(abs(pc2));

F_norm = F_relative(:);
F_norm = F_norm ./ max(abs(F_norm));

plot(t, pc1_n, 'k-', 'LineWidth',1.3, 'DisplayName','PC1');
hold on;
plot(t, pc2_n, 'r-', 'LineWidth',1.3, 'DisplayName','PC2');
plot(t, F_norm, 'k--', 'LineWidth',1.1, 'DisplayName','Flux (norm)');
yline(0, 'k:');
ylabel('Normalized units');
xlabel('Time (days)');
title(sprintf('PC1 & PC2 (%.1f%% / %.1f%% var)', explained(1), explained(2)));
legend('Location','best');
xlim([t(1) t(end)]);
grid on; box on;

% (c) Mode 1 (attenuation)
subplot(2,2,3);
imagesc(logD, depths, mode1);
set(gca,'YDir','reverse');
colormap(gca, turbo);
cb = colorbar;
cb.Label.String = 'EOF1 loading';
title(sprintf('Mode 1 (Attenuation, %.1f%%%% var)', explained(1)));
xlabel('D (mm)');
ylabel('Depth (m)');
set(gca,'XTick', log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});

% (d) Mode 2 (seesaw)
subplot(2,2,4);
imagesc(logD, depths, mode2);
set(gca,'YDir','reverse');
colormap(gca, turbo);
cb = colorbar;
cb.Label.String = 'EOF2 loading';
title(sprintf('Mode 2 (Seesaw, %.1f%%%% var)', explained(2)));
xlabel('D (mm)');
ylabel('Depth (m)');
set(gca,'XTick', log10([0.1 1 10]), 'XTickLabel',{'0.1','1','10'});

% --- FIGURE 2: 3D reconstructions ---------------------------------------
PSD_full_log = Y_log;
PSD_full_3d  = reshape(PSD_full_log, Nt, Nz, Ns);

Y1_log = (scores(:,1) * coeffs(:,1).') + Y_mean;
Y2_log = (scores(:,2) * coeffs(:,2).') + Y_mean;

Y1_3d = reshape(Y1_log, Nt, Nz, Ns);
Y2_3d = reshape(Y2_log, Nt, Nz, Ns);

% Colour limits from full model only (so modes compare cleanly)
vals = PSD_full_log(:);
vals = vals(isfinite(vals));
vals = sort(vals);
n    = numel(vals);
cmin = vals(max(1, round(0.02*n)));
cmax = vals(round(0.98*n));
clim_global = [cmin cmax];

figure('Name',[title_prefix ' 3D Reconstruction'], 'Color','w', ...
       'Position',[50 50 1100 950]);

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

% =========================================================================
% === HELPER: 3D "FENCE" PLOT =============================================
% =========================================================================
function plot_psd_slice_subplot(ax_handle, t_days, depths, D_mm, PSD_3d, ...
    title_str, c_limits)

axes(ax_handle); %#ok<LAXES>
hold on;

t_days = t_days(:);
depths = depths(:);
logD   = log10(D_mm(:));

% PSD_3d: [Nt x Nz x Ns] -> V: [Ns x Nt x Nz] for slice()
V = permute(PSD_3d, [3 1 2]);
[X, Y, Z] = meshgrid(t_days, logD, depths);

t_min = ceil(min(t_days));
t_max = floor(max(t_days));
xs    = t_min:1:t_max;        % one fence per day
xslice = fliplr(xs);          % earliest days in front

yslice = [];
zslice = [];

h = slice(X, Y, Z, V, xslice, yslice, zslice);
set(h, 'EdgeColor','none', 'FaceAlpha',0.95);

set(gca, 'ZDir','reverse');
colormap(gca, turbo);
caxis(c_limits);

cb = colorbar('Location','eastoutside');
cb.Label.String = 'log_{10} Part Vol (ppmV mm^{-1})';

xlabel('Time (days)');
ylabel('D (mm)');
zlabel('Depth (m)');
title(title_str);

set(gca,'YTick', log10([0.1 1 10]), 'YTickLabel',{'0.1','1','10'});

pbaspect([35 1 3]);
view([-60 5]);
axis tight;
box on; grid on;
end

% =========================================================================
% === INTERNAL ATLANTIC RUN + mu SENSITIVITY ==============================
% =========================================================================
function [out, cfg, eps_run, flux_struct] = run_simulation_internal()

if ~isfile('epsilon_daily.mat')
    error('Missing epsilon_daily.mat for internal EOF run.');
end

S = load('epsilon_daily.mat');
[t_days, eps_raw] = grab_data_simple(S);

mask    = t_days <= 60;
t_run   = t_days(mask);
eps_run = eps_raw(mask);

% --- base configuration (used for EOF, mu = 0.05) -----------------------
cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 200;
cfg.dz         = 10;

cfg.disagg_use_nonlinear = true;
cfg.disagg_kmax_a        = 0.95;
cfg.disagg_beta          = 1.0;
cfg.disagg_redistribute_p = 1.5;

cfg.attenuation_rate = 0.05;      % reference mu

cfg.use_NPP  = true;
cfg.NPP_rate = 2e-4;

% --- run reference simulation -------------------------------------------
sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_run, eps_run);
sim.run('tspan', t_run);

res = sim.result;
t   = res.time(:);
Y   = res.concentrations;

% compute export flux + relative for this run (for PCs plot)
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;

r_cm = sim.grid.getFractalRadii();
r_v  = sim.grid.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, sim.grid.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;
V_m3    = (4/3) * pi * r_v.^3 * 1e-6;

export_depth_idx = Nz;
Q_export = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);
F_total  = sum(Q_export .* (V_m3.' .* W_m_d.'), 2);

eps_floor = 1e-10;
idx_ss = find(eps_run < eps_floor);
if numel(idx_ss) < 5
    idx_ss = 1:min(10, numel(t));
end
F_ss = mean(F_total(idx_ss));
if F_ss == 0, F_ss = 1; end
F_relative = F_total / F_ss;

flux_struct.F_total    = F_total;
flux_struct.F_relative = F_relative;
flux_struct.F_ss       = F_ss;

% package outputs for EOF
out.t_days         = t;
out.concentrations = Y;
out.D_um           = 2 * r_cm * 1e4;

% --- OPTIONAL: attenuation sensitivity (multiple mu) --------------------
mu_list = [0.01 0.03 0.05 0.10];

figure('Name','Attenuation sensitivity (relative export flux)', ...
       'Color','w', 'Position',[200 200 700 450]);
hold on;

colors = lines(numel(mu_list));

for im = 1:numel(mu_list)
    mu_cfg = cfg;
    mu_cfg.attenuation_rate = mu_list(im);

    sim_mu = CoagulationSimulation(mu_cfg);
    sim_mu.setEpsilonTimeSeries(t_run, eps_run);
    sim_mu.run('tspan', t_run);

    Y_mu = sim_mu.result.concentrations;

    Q_export_mu = Y_mu(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);
    F_total_mu  = sum(Q_export_mu .* (V_m3.' .* W_m_d.'), 2);

    % normalize each series by its own steady state for clarity
    idx_ss_mu = idx_ss;
    F_ss_mu   = mean(F_total_mu(idx_ss_mu));
    if F_ss_mu == 0, F_ss_mu = 1; end
    F_rel_mu  = F_total_mu / F_ss_mu;

    plot(t_run, F_rel_mu, '-', 'Color', colors(im,:), ...
        'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\mu = %.2f d^{-1}', mu_list(im)));
end

xlabel('Time (days)');
ylabel('Relative export flux (F/F_{ss})');
title('Sensitivity of export flux to attenuation rate \mu');
grid on; box on;
legend('Location','best');

end

% =========================================================================
% === PLACEHOLDER FOR 0-D EXPERIMENT (NOT USED HERE) ======================
% =========================================================================
function onefile_obs_eps_allin
fprintf('Skipping 0-D sanity check in this version.\n');
end

% =========================================================================
% === SIMPLE LOADER FOR epsilon FILES ====================================
% =========================================================================
function [t, e] = grab_data_simple(S)
% Try to find time and epsilon vectors in a loosely structured .mat

names = fieldnames(S);
if numel(names) == 1 && isstruct(S.(names{1}))
    S = S.(names{1});
    names = fieldnames(S);
end

t = [];
e = [];

for k = 1:numel(names)
    val = S.(names{k});
    nm  = lower(names{k});

    if isnumeric(val)
        if (contains(nm,'time') || contains(nm,'day'))
            t = val;
        end
        if (contains(nm,'eps') || contains(nm,'dissip'))
            e = val;
        end
    end
end

if isempty(t) || isempty(e)
    error('Could not find time/epsilon arrays in epsilon file.');
end

t = t(:);
e = e(:);

% make time start at 0
t = t - t(1);

n = min(numel(t), numel(e));
t = t(1:n);
e = e(1:n);

% avoid zeros in epsilon
e = max(e, 1e-12);
end