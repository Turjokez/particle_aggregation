function plot_full_figure(res, od, gridObj, cfg, outdir, tag)
% plot_full_figure — compact 2×2 overview for one simulation
%
% Usage:
%   plot_full_figure(res, od, sim.grid, sim.config, outdir, 'baseline')

if nargin < 6, tag = 'run'; end
if ~exist(outdir,'dir'), mkdir(outdir); end

% ---------- time + epsilon(t)
t = res.time(:);

if isfield(res,'epsilon_used') && ~isempty(res.epsilon_used)
    eps_used = max(res.epsilon_used(:), 0);
else
    eps_used = zeros(size(t));
    for i = 1:numel(t)
        switch lower(cfg.epsilon_profile)
            case 'observed'
                if ~isempty(cfg.epsilon_time) && ~isempty(cfg.epsilon_series)
                    eps_used(i) = interp1(cfg.epsilon_time(:), cfg.epsilon_series(:), ...
                                           t(i), 'linear', 'extrap');
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

% ---------- bin diameters (robust)
if isfield(od,'D_cm')
    D_cm = od.D_cm(:);
elseif isfield(od,'diam_cm')
    D_cm = od.diam_cm(:);
elseif ~isempty(gridObj) && (ismethod(gridObj,'getFractalRadii') ...
        || any(strcmp(methods(gridObj),'getFractalRadii')))
    D_cm = 2 * gridObj.getFractalRadii();
else
    error('plot_full_figure:diameters','Could not determine bin diameters.');
end
D_cm = D_cm(:);
D_um = D_cm * 1e4;

% ---------- spectra matrix (t × sections)
N = [];
for nm = {'spectra','N','concentrations'}
    if isfield(od, nm{1}), N = od.(nm{1}); break; end
end
if isempty(N), N = res.concentrations; end
N = max(N, 0);
[nT, nS] = size(N); %#ok<NASGU>

% ---------- export flux from settling loss (with safe fallback)
flux_label = 'Settling loss (export) [model units d^{-1}]';
flux = [];
if isfield(res,'diagnostics') && isfield(res.diagnostics,'total_losses') ...
        && isfield(res.diagnostics.total_losses,'settl')
    flux = res.diagnostics.total_losses.settl(:);
end
if isempty(flux)
    % fallback: volume proxy (Σ N D^3). Label clearly as proxy.
    flux = sum(N .* reshape(D_cm(:)'.^3,1,[]), 2);
    flux_label = 'Volume proxy (\Sigma N D^3) [a.u.]';
end

% ---------- figure
fh = figure('Color','w','Position',[120 120 1100 700]);

% (A) ε(t) with threshold
ax1 = subplot(2,2,1);
plot(ax1, t, eps_used, '-', 'LineWidth',1.6); hold(ax1,'on');
yline(ax1, eps_ref, '--', 'LineWidth',1.0);
grid(ax1,'on'); set(ax1,'YScale','log');
xlabel(ax1,'Time (days)'); ylabel(ax1,'\epsilon (W kg^{-1})');
title(ax1, sprintf('\\epsilon(t) &  threshold  (%.1e)', eps_ref));
legend(ax1, {'\epsilon(t)','\epsilon_{ref}'}, 'Location','best');

% optional shaded storms if present
if isprop(cfg,'storm_windows') && ~isempty(cfg.storm_windows)
    yl = ylim(ax1);
    for r = 1:size(cfg.storm_windows,1)
        x = cfg.storm_windows(r,:);
        patch(ax1, [x(1) x(2) x(2) x(1)], [yl(1) yl(1) yl(2) yl(2)], ...
              [0.85 0.9 1], 'EdgeColor','none','FaceAlpha',0.35);
    end
    uistack(findobj(ax1,'Type','line'),'top');
end

% (B) Export (true settling loss if available)
ax2 = subplot(2,2,2);
plot(ax2, t, flux, '-', 'LineWidth',1.6);
grid(ax2,'on');
xlabel(ax2,'Time (days)'); ylabel(ax2, flux_label);
title(ax2,'Export time series');

% (C) PSD at start vs end
ax3 = subplot(2,2,3);
loglog(ax3, D_um, max(N(1,:),1e-20), '-', 'LineWidth',1.6); hold(ax3,'on');
loglog(ax3, D_um, max(N(end,:),1e-20), '-', 'LineWidth',1.6);
grid(ax3,'on');
xlabel(ax3,'Diameter (\mum)'); ylabel(ax3,'# cm^{-3}');
legend(ax3, {sprintf('t=%.1f d', t(1)), sprintf('t=%.1f d', t(end))}, 'Location','southwest');
title(ax3,'PSD at start vs end');

% (D) Size–time heatmap of log10(N)
ax4 = subplot(2,2,4);
imagesc(ax4, t, log10(D_um(:)), log10(max(N,1e-20))'); axis(ax4,'xy');
xlabel(ax4,'Time (days)'); ylabel(ax4,'log_{10} D (\mum)');
title(ax4,'log_{10} PSD (time × size)');
cb = colorbar(ax4); cb.Label.String = 'log_{10} # cm^{-3}';
try
    colormap(ax4, turbo(256));
catch
    colormap(ax4, parula(256));
end

sgtitle(sprintf('%s — %s', strrep(tag,'_','\_'), datestr(now,'yyyy-mm-dd HH:MM')));

% remove axes toolbars in exports
try
    axs = findall(fh,'Type','axes'); for k=1:numel(axs), axs(k).Toolbar = []; end
end

% ---------- save
pngf = fullfile(outdir, sprintf('full_figure_%s.png', tag));
figf = fullfile(outdir, sprintf('full_figure_%s.fig', tag));
exportgraphics(fh, pngf, 'Resolution', 200);
savefig(fh, figf);
fprintf('✔ Full figure saved:\n  %s\n  %s\n', pngf, figf);
end