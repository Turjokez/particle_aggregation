function mu_NPP_diagnostics
% Simple mu–NPP diagnostics to understand why PSD goes "all red" at the end.

outdir = 'fig_muNPP_tests';
if ~exist(outdir,'dir'); mkdir(outdir); end

if ~isfile('epsilon_daily.mat')
    error('epsilon_daily.mat not found in current folder.');
end
S = load('epsilon_daily.mat');
[t_days, eps_raw] = grab_data_simple(S);

cases(1).label    = 'strongMu_noNPP';
cases(1).title    = 'Strong attenuation, no NPP';
cases(1).mu       = 0.30;
cases(1).use_NPP  = false;
cases(1).NPP_rate = 0.0;
cases(1).t_final  = 15;

cases(2).label    = 'weakMu_highNPP';
cases(2).title    = 'Weak attenuation, high NPP';
cases(2).mu       = 0.01;
cases(2).use_NPP  = true;
cases(2).NPP_rate = 2e-4;
cases(2).t_final  = 15;

cases(3).label    = 'atlanticLike_tuned';
cases(3).title    = 'Atlantic-like tuned (mu=0.10, NPP=5e-5)';
cases(3).mu       = 0.10;
cases(3).use_NPP  = true;
cases(3).NPP_rate = 5e-5;
cases(3).t_final  = 30;

for ic = 1:numel(cases)
    C = cases(ic);
    fprintf('\n=== CASE %d: %s ===\n', ic, C.title);

    mask    = t_days <= C.t_final;
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

    cfg.attenuation_rate = C.mu;
    cfg.use_NPP          = C.use_NPP;
    cfg.NPP_rate         = C.NPP_rate;

    cfg.t_init  = t_run(1);
    cfg.t_final = t_run(end);
    cfg.delta_t = t_run(2) - t_run(1);

    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    result = sim.run('tspan', t_run);

    [f1, f2] = make_case_plots(result, sim.grid, cfg, eps_run, C.title);

    fname1 = fullfile(outdir, sprintf('%02d_%s_flux_inventory.png', ic, C.label));
    fname2 = fullfile(outdir, sprintf('%02d_%s_PSD3D.png',          ic, C.label));

    exportgraphics(f1, fname1, 'Resolution',300);
    exportgraphics(f2, fname2, 'Resolution',300);

    close(f1); close(f2);
end

fprintf('\nAll diagnostic figures saved in folder: %s\n', outdir);
end

%% ------------------------------------------------------------------------
function [f1, f2] = make_case_plots(result, grid_obj, cfg, eps_run, cfg_label)

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

export_depth_idx = Nz;
Q_export = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);
F_total = sum(Q_export .* V_m3.' .* W_m_d.', 2);

F_ref_idx = 1:min(10,Nt);
F_ref     = mean(F_total(F_ref_idx));
if F_ref == 0, F_ref = 1; end
F_rel = F_total / F_ref;

inventory_simple = sum(Y, 2);

Y_3d = reshape(Y, Nt, Nz, Ns);
inventory_vol = zeros(Nt,1);
for it = 1:Nt
    slice_it = squeeze(Y_3d(it,:,:));
    inventory_vol(it) = sum(slice_it * V_m3(:));
end

inv_simple_rel = inventory_simple / max(inventory_simple(1),1e-30);
inv_vol_rel    = inventory_vol    / max(inventory_vol(1),1e-30);

f1 = figure('Name',[cfg_label ' – Flux & Inventory'], ...
            'Color','w','Position',[100 100 800 900]);

subplot(3,1,1);
semilogy(t, max(eps_run,1e-12),'k-','LineWidth',1.2); hold on;
plot(t([1 end]),[1e-6 1e-6],'r--','LineWidth',1);
ylabel('\epsilon (W kg^{-1})');
title(['Turbulence Forcing - ' cfg_label]);
grid on; axis tight;

subplot(3,1,2);
yyaxis left;
plot(t, F_total,'b-','LineWidth',2);
ylabel('Export Flux (m^3 d^{-1})');
yyaxis right;
plot(t, F_rel,'Color',[0.85 0.4 0],'LineWidth',1.5);
ylabel('Relative Flux (F/F_{ref})');
title('Export Flux at Bottom Layer');
xlabel('Time (days)');
grid on; axis tight;

subplot(3,1,3);
plot(t, inv_simple_rel,'r-','LineWidth',2,'DisplayName','Number inventory'); hold on;
plot(t, inv_vol_rel,'b-','LineWidth',2,'DisplayName','Volume inventory');
yline(1,'k--','LineWidth',1);
xlabel('Time (days)');
ylabel('Normalized Inventory');
title('Total Particle Inventory');
legend('Location','best');
grid on; axis tight;

D_um  = 2 * r_cm(:) * 1e4;
D_mm  = D_um / 1000;
depths = (1:Nz) * cfg.dz;

Y_log  = log10(max(Y, 1e-20));
PSD_3d = reshape(Y_log, Nt, Nz, Ns);

f2 = figure('Name',[cfg_label ' – 3D PSD Full Field'], ...
            'Color','w','Position',[100 100 1000 800]);
ax = axes;
plot_psd_slice_simple(ax, t, depths, D_mm, PSD_3d, ...
    'Full Model (log_{10} Part Vol)');
end

%% ------------------------------------------------------------------------
function plot_psd_slice_simple(ax_handle, t_days, depths, D_mm, PSD_3d, title_str)

axes(ax_handle); hold on;

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

vals = PSD_3d(:);
vals = vals(isfinite(vals));
if ~isempty(vals)
    n  = numel(vals);
    lo = vals(max(1, round(0.02*n)));
    hi = vals(round(0.98*n));
    if hi <= lo, hi = lo + 1; end
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

%% ------------------------------------------------------------------------
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