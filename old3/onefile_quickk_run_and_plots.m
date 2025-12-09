function onefile_quickk_run_and_plots()
% ONEFILE_QUICKK_RUN_AND_PLOTS (robust)
% - Loads observed turbulence (epsilon) from your MAT file (case-insensitive).
% - If no time field is present, synthesizes t_days = 0:(NT-1).
% - Handles ε(t) or ε(z,t) with optional MLD-mean.
% - Runs the coagulation model with epsilon(t).
% - Saves a compact overview figure panel.

    % -------- Paths / I-O --------
    MATFILE = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code/epsilon_daily.mat';
    if ~exist(MATFILE,'file')
        error('MAT file not found:\n  %s\nEdit MATFILE in this script.', MATFILE);
    end
    ts = datestr(now,'yyyymmdd_HHMMSS');
    outdir = fullfile('/Users/turjo/Desktop/run_obs_eps_master', ts);
    if ~exist(outdir,'dir'), mkdir(outdir); end
    tag = 'observed_eps';

    fprintf('=== onefile_quickk_run_and_plots ===\n');
    fprintf('MAT file: %s\n', MATFILE);
    fprintf('Out dir : %s\n', outdir);

    % -------- Load epsilon(t) --------
    S0 = load(MATFILE);
    S  = unwrap_singleton_structs(S0);
    [t_days, eps_series] = grab_time_eps(S);
    fprintf('Time span: %.2f to %.2f days (%d samples)\n', t_days(1), t_days(end), numel(t_days));
    fprintf('ε range  : [%.2e, %.2e]\n', min(eps_series), max(eps_series));

    % -------- Configure & run --------
    sim_days = max(1, floor(t_days(end)));
    cfg = SimulationConfig( ...
        'epsilon_profile','observed', ...
        'epsilon_time', t_days, ...
        'epsilon_series', eps_series, ...
        'epsilon_ref', 1e-6, ...
        'p_alpha', 0.2, ...                % mild α(ε) scaling
        'disagg_use_nonlinear', true, ...  % external breakup on
        't_init', 0, ...
        't_final', sim_days, ...
        'delta_t', 1 ...
    );

    sim = CoagulationSimulation(cfg);
    res = sim.run();
    sim.generateOutputs(true);

    % -------- Figure panel --------
    od = res.output_data;
    plot_full_figure(res, od, sim.grid, sim.config, outdir, tag);
end

% ======================================================================
% Helpers
% ======================================================================
function S = unwrap_singleton_structs(S)
    % If S is a struct with exactly one field that is itself a struct,
    % unwrap it (do this repeatedly).
    while isstruct(S)
        f = fieldnames(S);
        if numel(f)==1 && isstruct(S.(f{1}))
            S = S.(f{1});
        else
            break;
        end
    end
end

function [t_days, eps_series] = grab_time_eps(S)
% Robust detection of time and epsilon from struct S (case-insensitive).
% Returns:
%   t_days     - column vector in days from 0
%   eps_series - column vector ε(t)

    % ---- Make a case-insensitive field index ----
    fn = fieldnames(S);
    fn_low = lower(fn);
    map = containers.Map(fn_low, fn); % lower->original

    % ---- Helper to fetch first matching field (case-insensitive) ----
    function name = first_ci(cands)
        name = '';
        for i=1:numel(cands)
            key = lower(cands{i});
            if isKey(map, key), name = map(key); return; end
        end
    end

    % ---- epsilon field ----
    e_name = first_ci({'eps','epsilon','epsilon_daily','eps_ml','eps_series','E'});
    if isempty(e_name)
        % fallback: any 1-D/2-D numeric likely to be epsilon
        nameGuess = guess_numeric_timevarying_field(S, map);
        if isempty(nameGuess)
            error('Could not find epsilon field.');
        end
        e_name = nameGuess;
    end
    E = S.(e_name);

    % ---- time field (try many) ----
    t_name = first_ci({'mtime','time','t_days','days','t','datenum','datetime','date_time','dateTime','serial_time'});
    t_days = [];
    if ~isempty(t_name)
        t_raw = S.(t_name);
        t_days = normalize_time_to_days(t_raw);
        fprintf('grab_time_eps: using time field "%s" (N=%d)\n', t_name, numel(t_days));
    end

    % ---- depth/MLD (optional) ----
    z_name   = first_ci({'z','depth','Z'});
    mld_name = first_ci({'mld','MLD','mixedlayerdepth','mixed_layer_depth'});

    % ---- build eps_series ----
    if isvector(E)
        % 1-D epsilon
        eps_series = E(:);
        if isempty(t_days)
            t_days = (0:numel(eps_series)-1).'; %#ok<*NBRAK>
            fprintf('grab_time_eps: synthesized time 0..%d (N=%d) from ε length.\n', numel(eps_series)-1, numel(eps_series));
        elseif numel(eps_series) ~= numel(t_days)
            % resample ε to time length
            eps_series = interp1(linspace(0,1,numel(eps_series)), eps_series, linspace(0,1,numel(t_days)),'linear','extrap').';
            fprintf('grab_time_eps: resampled ε to match time length.\n');
        end

    elseif isnumeric(E) && ismatrix(E)
        % 2-D epsilon (assume ε(z,t) or ε(t,z))
        % Decide time dimension: choose the dimension that matches time length, else the larger one
        [nz, nt] = size(E);
        tdim = 2; % default assume columns are time
        if ~isempty(t_days)
            if numel(t_days)==nz, tdim = 1; end
            if numel(t_days)==nt, tdim = 2; end
        else
            % no time: pick longer dimension as time
            if nt >= nz, tdim = 2; else, tdim = 1; end
            Nguess = (tdim==2) * nt + (tdim==1) * nz;
            t_days = (0:Nguess-1).';
            fprintf('grab_time_eps: synthesized time 0..%d (N=%d) from ε matrix.\n', Nguess-1, Nguess);
        end

        if tdim == 2
            % E(:,t)
            if ~isempty(z_name) && ~isempty(mld_name)
                z   = S.(z_name);   z = z(:);
                mld = S.(mld_name); mld = mld(:);
                mld = match_length(mld, numel(t_days));
                eps_series = mld_mean(E, z, mld);
                fprintf('grab_time_eps: using epsilon field "%s" (size [%d %d])\n', e_name, nz, nt);
            else
                k = min(3, nz);
                eps_series = squeeze(mean(E(1:k, :), 1)).';
                fprintf('grab_time_eps: depth info missing — mean of top-%d levels.\n', k);
            end
        else
            % E(t,:)
            if ~isempty(z_name) && ~isempty(mld_name)
                z   = S.(z_name);   z = z(:);
                mld = S.(mld_name); mld = match_length(mld, numel(t_days));
                eps_series = mld_mean(E.', z, mld); % transpose so columns are time
                fprintf('grab_time_eps: using epsilon field "%s" (size [%d %d], transposed)\n', e_name, nz, nt);
            else
                k = min(3, nt);
                eps_series = squeeze(mean(E(:, 1:k), 2));
                eps_series = eps_series(:);
                fprintf('grab_time_eps: depth info missing — mean of first-%d columns (assumed top levels).\n', k);
            end
        end
    else
        error('Unsupported epsilon array type/shape.');
    end

    % Smooth & sanitize
    if numel(eps_series) >= 3
        eps_series = movmean(eps_series, 3, 'Endpoints','shrink');
    end
    eps_series = max(eps_series, 0);
    eps_series = max(eps_series, 1e-12);  % floor to avoid α blow-ups
    eps_series = eps_series(:);
    t_days     = t_days(:);
end

function t_days = normalize_time_to_days(t_raw)
    if isdatetime(t_raw)
        t_days = days(t_raw(:) - t_raw(1));
        return;
    end
    t_raw = t_raw(:);
    % Heuristics: datenum ~ 7e5; epoch sec ~ 1e9; epoch days ~ 2e4
    tr = double(t_raw);
    if max(tr) > 1e8                  % likely epoch seconds
        t_days = (tr - tr(1))/86400;
    elseif max(tr) > 1e5              % likely datenum
        t_days = tr - tr(1);
    else                               % likely already in days
        t_days = tr - tr(1);
    end
end

function v = match_length(v, N)
    v = v(:);
    if numel(v) == N, return; end
    v = interp1(linspace(0,1,numel(v)), v, linspace(0,1,N), 'linear', 'extrap').';
end

function eps_ml = mld_mean(Ezt, z, mld)
    % Ezt is ε(z,t) with columns as time
    nz = size(Ezt,1);
    eps_ml = zeros(numel(mld),1);
    for j = 1:numel(mld)
        mask = z <= mld(j);
        if ~any(mask)
            k = min(3, nz);
            eps_ml(j) = mean(Ezt(1:k, j), 'omitnan');
        else
            eps_ml(j) = mean(Ezt(mask, j), 'omitnan');
        end
    end
end

function nameGuess = guess_numeric_timevarying_field(S, map)
    % Try to find a numeric field with 1D or 2D size (time-varying)
    nameGuess = '';
    keys = map.keys;
    for i=1:numel(keys)
        fname = map(keys{i});
        v = S.(fname);
        if isnumeric(v) && (isvector(v) || ismatrix(v))
            % Skip obviously non-ε things: lat/lon with small ranges
            if isvector(v)
                if numel(v) >= 12 && range(double(v)) > 0
                    nameGuess = fname; return;
                end
            else
                if all(size(v) >= 12)
                    nameGuess = fname; return;
                end
            end
        end
    end
end

% ======================================================================
% Figure panel (arg name gridObj avoids clashing with MATLAB grid())
% ======================================================================
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