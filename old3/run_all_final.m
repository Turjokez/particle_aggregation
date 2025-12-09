%% onefile_quick_run_and_plots.m
% Single-file driver + plotting + metrics (no extra .m files needed)

%% --- PATHS & INPUTS ------------------------------------------------------
project_root = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation';
addpath(genpath(fullfile(project_root,'src')));   % classes already in your repo

matfile = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code/epsilon_daily.mat';
outdir  = '/Users/turjo/Desktop/run_obs_eps_master/quick_onefile';
if ~exist(outdir,'dir'), mkdir(outdir); end

% Load observed epsilon series (replace with your loader if different)
% Expected fields inside MAT: mtime (datenum or days), eps (W/kg)
S = load(matfile); 
if isfield(S,'mtime'), t_days = S.mtime(:); else, error('Need S.mtime'); end
if isfield(S,'eps'),   eps_series = S.eps(:); else, error('Need S.eps'); end

% Align to model-day origin (t_init = 0 at first day)
t_days = t_days - t_days(1);

%% --- BASE CONFIG ---------------------------------------------------------
cfg = SimulationConfig( ...
    'epsilon_profile','observed', ...
    'epsilon_time',   t_days, ...
    'epsilon_series', eps_series, ...
    'epsilon_ref',    1e-6, ...      % 10^{-6} W/kg threshold line
    'p_alpha',        0.0,  ...      % start with no alpha(ε) scaling
    'disagg_use_nonlinear', true, ...
    't_init', 0, 't_final', 30, 'delta_t', 1);

%% --- RUN 1: baseline (p_alpha = 0) --------------------------------------
fprintf('\n=== RUN: baseline (p_\\alpha=%.2f) ===\n', cfg.p_alpha);
sim = CoagulationSimulation(cfg);
res = sim.run();
sim.generateOutputs(true);  % your existing full outputs

od = res.output_data;
plot_full_figure(res, od, sim.grid, sim.config, outdir, sprintf('baseline_pa%.2g',cfg.p_alpha));
summarize_run_metrics(res, od);

%% --- RUN 2 & 3: small sweep of p_alpha -----------------------------------
pAlphaList = [0.2, 0.4];
for pa = pAlphaList
    cfg.p_alpha = pa;
    fprintf('\n=== RUN: p_\\alpha=%.2f ===\n', pa);
    sim = CoagulationSimulation(cfg);
    res = sim.run();
    sim.generateOutputs(false);  % skip the heavy plotAll

    od = res.output_data;
    tag = sprintf('pa%.2g', pa);
    plot_full_figure(res, od, sim.grid, sim.config, outdir, tag);
    summarize_run_metrics(res, od);
end

fprintf('\nAll done. PNGs in: %s\n', outdir);

%% ========================================================================
%% =============== Local functions (stay in THIS file) ====================
%% ========================================================================

function plot_full_figure(res, od, grid, cfg, outdir, tag)
% plot_full_figure — compact overview panel for one simulation (local)

if nargin < 6, tag = 'run'; end
if ~exist(outdir,'dir'), mkdir(outdir); end

% --- time & epsilon
t = res.time(:);
if isfield(res,'epsilon_used') && ~isempty(res.epsilon_used)
    eps_used = max(res.epsilon_used(:),0);
else
    eps_used = zeros(size(t));
    for i = 1:numel(t)
        switch lower(cfg.epsilon_profile)
            case 'observed'
                if ~isempty(cfg.epsilon_time) && ~isempty(cfg.epsilon_series)
                    eps_used(i) = interp1(cfg.epsilon_time(:), cfg.epsilon_series(:), t(i), 'linear', 'extrap');
                else
                    eps_used(i) = cfg.epsilon;
                end
            case 'sine'
                eps_used(i) = cfg.epsilon_mean + cfg.epsilon_amp * ...
                              sin(2*pi*t(i)/cfg.epsilon_period + cfg.epsilon_phase);
            otherwise
                eps_used(i) = cfg.epsilon;
        end
    end
end
eps_ref = cfg.epsilon_ref;

% --- diameters
D_cm = [];
if isfield(od,'D_cm')
    D_cm = od.D_cm(:);
elseif isfield(od,'diam_cm')
    D_cm = od.diam_cm(:);
elseif ismethod(grid,'getFractalRadii')
    r_cm = grid.getFractalRadii();   % <-- call first
    D_cm = 2 * r_cm(:);              % then reshape/index
else
    error('Could not determine bin diameters.');
end
D_um = D_cm * 1e4;

% --- spectra (time × size)
N = [];
for k = {'spectra','N','concentrations'}
    if isfield(od,k{1}), N = od.(k{1}); break; end
end
if isempty(N), N = res.concentrations; end
N = max(N,0);
[nT,~] = size(N);

% --- flux (if absent, proxy)
flux = [];
for fn = {'export_flux','flux','F'}
    if isfield(od, fn{1}), flux = od.(fn{1}); break; end
end
if isempty(flux)
    flux = sum(N .* (D_cm(:)'.^3), 2); % crude proxy ~ volume
end
flux = flux(:);

% --- figure
fh = figure('Color','w','Position',[120 120 1100 700]);

% (A) ε(t)
ax1 = subplot(2,2,1);
plot(ax1, t, eps_used, '-', 'LineWidth',1.6); hold(ax1,'on');
if ~isempty(eps_ref) && eps_ref>0, yline(ax1, eps_ref, '--', 'LineWidth',1.0); end
grid(ax1,'on'); ax1.YScale = 'log';
xlabel(ax1,'Time (days)'); ylabel(ax1,'\epsilon (W kg^{-1})');
ttl = '\epsilon(t)';
if ~isempty(eps_ref) && eps_ref>0, ttl = sprintf('%s & \\epsilon_{ref}=%.1e',ttl,eps_ref); end
title(ax1, ttl);
legend(ax1, {'\epsilon(t)','\epsilon_{ref}'}, 'Location','best');

% (B) Export flux
ax2 = subplot(2,2,2);
plot(ax2, t, flux, '-', 'LineWidth',1.6);
grid(ax2,'on');
xlabel(ax2,'Time (days)'); ylabel(ax2,'Export flux (a.u.)');
title(ax2,'Export flux');

% (C) PSD start vs end
ax3 = subplot(2,2,3);
loglog(ax3, D_um, max(N(1,:),1e-20), '-', 'LineWidth',1.6); hold(ax3,'on');
loglog(ax3, D_um, max(N(nT,:),1e-20), '-', 'LineWidth',1.6);
grid(ax3,'on'); xlabel(ax3,'Diameter (\mum)'); ylabel(ax3,'# cm^{-3}');
legend(ax3, {sprintf('t=%.1f d', t(1)), sprintf('t=%.1f d', t(end))}, 'Location','southwest');
title(ax3,'PSD at start vs end');

% (D) Heatmap log10 N
ax4 = subplot(2,2,4);
imagesc(ax4, t, log10(D_um(:)), log10(max(N,1e-20))'); axis(ax4,'xy');
xlabel(ax4,'Time (days)'); ylabel(ax4,'log_{10} D (\mum)');
title(ax4,'log_{10} PSD (time × size)');
cb = colorbar(ax4); cb.Label.String = 'log_{10} # cm^{-3}';
colormap(ax4, turbo(256));

sgtitle(sprintf('%s — %s', tag, datestr(now,'yyyy-mm-dd HH:MM')));

% save
pngf = fullfile(outdir, sprintf('full_figure_%s.png', tag));
figf = fullfile(outdir, sprintf('full_figure_%s.fig', tag));
exportgraphics(fh, pngf, 'Resolution', 200);
savefig(fh, figf);
fprintf('✔ Full figure saved:\n  %s\n  %s\n', pngf, figf);
end

function summarize_run_metrics(res, od)
% summarize_run_metrics — quick text summary (local)

eps_used = []; 
if isfield(res,'epsilon_used'), eps_used = res.epsilon_used(:); end
if ~isempty(eps_used)
    fprintf('ε used: min=%.2e, median=%.2e, max=%.2e\n', ...
        min(eps_used), median(eps_used), max(eps_used));
end

% crude τ50 proxy from cumulative flux if present
if isfield(od,'export_flux')
    F = od.export_flux(:);
    CF = cumsum(F); 
    if CF(end) > 0
        idx = find(CF >= 0.5*CF(end), 1, 'first');
        fprintf('τ50 (time to half export): t=%.2f d\n', res.time(idx));
    end
end

% simple ε–flux correlation if both exist
flux = [];
for fn = {'export_flux','flux','F'}
    if isfield(od, fn{1}), flux = od.(fn{1}); break; end
end
if ~isempty(eps_used) && ~isempty(flux)
    N = min(numel(eps_used), numel(flux));
    r = corr(eps_used(1:N), flux(1:N), 'rows','complete', 'type','Spearman');
    fprintf('Spearman corr(ε, flux) = %.3f\n', r);
end
end