function master_simulation_runner2
% MASTER_SIMULATION_RUNNER2
%
% Diagnostic script to test attenuation vs NPP and check whether
% total particle inventory keeps growing ("all red at the end").
%
% Uses:
%   - SimulationConfig.m
%   - CoagulationSimulation.m
%   - SettlingVelocityService.m
%
% This script does NOT change the main project code; it just
% creates a simpler playground to tune attenuation / NPP.

clc;
fprintf('======================================================\n');
fprintf('   PARTICLE AGGREGATION MODEL - DIAGNOSTIC EDITION\n');
fprintf('======================================================\n');
fprintf('Choose a diagnostic configuration:\n');
fprintf('  [1] Short run, STRONG attenuation (mu = 0.3 d^{-1}), no NPP\n');
fprintf('  [2] Short run, WEAK attenuation (mu = 0.01 d^{-1}), with NPP\n');
fprintf('  [3] Atlantic-like, TUNED (mu = 0.10, NPP = 5e-05)\n');
fprintf('======================================================\n');

choice = input('Enter number [1-3]: ');

switch choice
    case 1
        cfg_label        = 'Strong attenuation, no NPP';
        params.mu        = 0.30;     % strong loss
        params.use_NPP   = false;
        params.NPP_rate  = 0.0;
        params.t_final   = 15;       % shorter run
    case 2
        cfg_label        = 'Weak attenuation, with NPP';
        params.mu        = 0.01;     % weak loss
        params.use_NPP   = true;
        params.NPP_rate  = 2e-4;
        params.t_final   = 15;       % shorter run
    case 3
        cfg_label        = 'Atlantic-like (mu = 0.10, NPP = 5e-05)';
        params.mu        = 0.10;     % *** TUNED ATTENUATION ***
        params.use_NPP   = true;
        params.NPP_rate  = 5e-5;     % *** TUNED NPP ***
        params.t_final   = 30;       % full Atlantic eddy window
    otherwise
        fprintf('Invalid selection.\n');
        return;
end

fprintf('\n=== RUNNING: %s ===\n', cfg_label);

% -------------------------------------------------------------------------
% 1. Load epsilon forcing (same Atlantic forcing as main script)
% -------------------------------------------------------------------------
if ~isfile('epsilon_daily.mat')
    error('Missing epsilon_daily.mat in current folder.');
end

S = load('epsilon_daily.mat');
[t_days, eps_raw] = grab_data_simple_local(S);

% Restrict to chosen time window
mask    = t_days <= params.t_final;
t_run   = t_days(mask);
eps_run = eps_raw(mask);

fprintf('Simulation window: %.1f days (N = %d time steps)\n', ...
    t_run(end) - t_run(1), numel(t_run));

% -------------------------------------------------------------------------
% 2. Build SimulationConfig
% -------------------------------------------------------------------------
cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max      = 200;
cfg.dz         = 10;

% Nonlinear disaggregation (same as main setup)
cfg.disagg_use_nonlinear  = true;
cfg.disagg_kmax_a         = 0.95;
cfg.disagg_beta           = 1.0;
cfg.disagg_redistribute_p = 1.5;

% Attenuation + NPP knobs
cfg.attenuation_rate = params.mu;
cfg.use_NPP          = params.use_NPP;
cfg.NPP_rate         = params.NPP_rate;

% Time settings
cfg.t_init  = t_run(1);
cfg.t_final = t_run(end);
cfg.delta_t = t_run(2) - t_run(1);   % assume uniform spacing

% -------------------------------------------------------------------------
% 3. Run the model
% -------------------------------------------------------------------------
sim = CoagulationSimulation(cfg);
sim.setEpsilonTimeSeries(t_run, eps_run);

result = sim.run('tspan', t_run);

% -------------------------------------------------------------------------
% 4. Diagnostics: flux + inventory + 3D PSD
% -------------------------------------------------------------------------
diagnostic_plots(result, sim.grid, cfg, eps_run, cfg_label);

end

%% ========================================================================
%  SIMPLE EPSILON LOADER (LOCAL COPY)
% ========================================================================
function [t, e] = grab_data_simple_local(S)

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
e = max(e, 1e-12);   % floor to avoid zeros
end

%% ========================================================================
%  DIAGNOSTIC PLOTS:
%   1) epsilon(t)
%   2) export flux (absolute + relative)
%   3) total inventory
%   4) 3D PSD fence plot
% ========================================================================
function diagnostic_plots(result, grid_obj, cfg, eps_run, cfg_label)

t  = result.time(:);            % Nt x 1
Y  = result.concentrations;     % Nt x (Nz*Ns)
Nt = numel(t);

Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;

% -------------------------------------------------------------------------
% A) Export flux at bottom
% -------------------------------------------------------------------------
r_cm    = grid_obj.getFractalRadii();
r_v     = grid_obj.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_obj.setcon);

% Convert
W_m_d = (ws_cm_s / 100) * 86400;      % cm/s -> m/day
V_m3  = (4/3) * pi * r_v.^3 * 1e-6;   % cm^3 -> m^3

export_depth_idx = Nz;
Q_export = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);   % Nt x Ns

F_total = sum(Q_export .* V_m3.' .* W_m_d.', 2);   % Nt x 1

% Relative to initial low-turbulence mean
F_ref_idx = 1:min(10, Nt);
F_ref     = mean(F_total(F_ref_idx));
if F_ref == 0, F_ref = 1; end
F_rel = F_total / F_ref;

% -------------------------------------------------------------------------
% B) Total inventory (simple and volume-weighted)
% -------------------------------------------------------------------------
% Simple inventory = sum of all state variables
inventory_simple = sum(Y, 2);           % Nt x 1

% Volume-weighted inventory (sum over depth and size of C * volume)
Y_3d = reshape(Y, Nt, Nz, Ns);          % [time, depth, size]
inventory_vol = zeros(Nt, 1);
for it = 1:Nt
    slice_it = squeeze(Y_3d(it, :, :));     % [Nz x Ns]
    % slice_it: [Nz x Ns], V_m3: [Ns x 1]
    inventory_vol(it) = sum(slice_it * V_m3(:));   % scalar
end

inv_simple_rel = inventory_simple / max(inventory_simple(1), 1e-30);
inv_vol_rel    = inventory_vol    / max(inventory_vol(1),    1e-30);

% -------------------------------------------------------------------------
% C) FIGURE 1: epsilon, flux, inventory
% -------------------------------------------------------------------------
figure('Name',[cfg_label ' - Flux & Inventory'], ...
       'Color','w','Position',[100 100 800 900]);

% (1) epsilon(t)
subplot(3,1,1);
semilogy(t, max(eps_run,1e-12), 'k-','LineWidth',1.2); hold on;
plot(t([1 end]), [1e-6 1e-6], 'r--','LineWidth',1);
ylabel('\epsilon (W kg^{-1})');
title(['Turbulence Forcing - ' cfg_label]);
grid on; axis tight;

% (2) export flux
subplot(3,1,2);
yyaxis left;
plot(t, F_total, 'b-','LineWidth',2);
ylabel('Export Flux (m^3 d^{-1})');
yyaxis right;
plot(t, F_rel, 'k--','LineWidth',1.5);
ylabel('Relative Flux (F / F_{ref})');
title('Export Flux at Bottom Layer');
xlabel('Time (days)');
grid on; axis tight;

% (3) inventory
subplot(3,1,3);
plot(t, inv_simple_rel, 'r-','LineWidth',2,'DisplayName','Number inventory');
hold on;
plot(t, inv_vol_rel, 'b-','LineWidth',2,'DisplayName','Volume inventory');
yline(1,'k--','LineWidth',1);
xlabel('Time (days)');
ylabel('Normalized Inventory');
title('Total Particle Inventory (diagnostic)');
legend('Location','best');
grid on; axis tight;

% -------------------------------------------------------------------------
% D) FIGURE 2: 3D PSD fence plot (full field only)
% -------------------------------------------------------------------------
D_um  = 2 * r_cm(:) * 1e4;     % ESD in microns
D_mm  = D_um / 1000;
depths = (1:Nz) * cfg.dz;

Y_log  = log10(max(Y, 1e-20));
PSD_3d = reshape(Y_log, Nt, Nz, Ns);

figure('Name',[cfg_label ' - 3D PSD Full Field'], ...
       'Color','w','Position',[100 100 1000 800]);

ax = axes; %#ok<LAXES>
plot_psd_slice_simple(ax, t, depths, D_mm, PSD_3d, ...
    'Full Model (log_{10} Part Vol)');

end

%% ========================================================================
%  SIMPLE 3D FENCE PLOT (no EOF, just full field)
% ========================================================================
function plot_psd_slice_simple(ax_handle, t_days, depths, D_mm, PSD_3d, title_str)

axes(ax_handle); %#ok<LAXES>
hold on;

t_days = t_days(:);
depths = depths(:);
logD   = log10(D_mm(:));

V = permute(PSD_3d, [3 1 2]);     % [size, time, depth]
[X,Y,Z] = meshgrid(t_days, logD, depths);

t_min = ceil(min(t_days));
t_max = floor(max(t_days));
xs    = t_min:1:t_max;            % one slice per day
xslice = fliplr(xs);
yslice = [];
zslice = [];

h = slice(X,Y,Z,V,xslice,yslice,zslice);
set(h,'EdgeColor','none','FaceAlpha',0.95);

set(gca,'ZDir','reverse');
colormap(gca, turbo);

vals = PSD_3d(:);
vals = vals(isfinite(vals));
if ~isempty(vals)
    n  = numel(vals);
    lo = vals(max(1, round(0.02*n)));
    hi = vals(round(0.98*n));
    if hi <= lo
        hi = lo + 1;
    end
    caxis([lo hi]);
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