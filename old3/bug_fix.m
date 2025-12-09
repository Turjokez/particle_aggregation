function onefile_obs_eps_allin()
% ========================================================================
% onefile_obs_eps_allin.m  —  ε(t) → PSD → export (robust, self-contained)
%
% What this script does
%   1) Loads epsilon time-series from a .mat (robust field/shape handling)
%   2) Resamples ε(t) to coarse/fine Δt for robustness checks
%   3) Runs CoagulationSimulation (even if no ε hook is available)
%   4) Extracts absolute PSD (no per-time normalization)
%   5) Computes mass flux from PSD × w(D) (units g cm^-2 d^-1)
%   6) Makes overview, QA, and Δt robustness figures
%
% Notes
%   • If your class lacks a public epsilon hook, the script still runs and
%     computes export from PSD×w(D). The log will say: ε forcing hook used: 0
%   • No LaTeX errors: all labels use 'tex' interpreter and colorbar via ylabel(cb,...)
%   • Change velocity law in the PARAMS section below if desired.
% ========================================================================

fprintf('=== onefile_obs_eps_allin ===\n');

% ------------------------------------------------------------------------
% USER INPUTS
% ------------------------------------------------------------------------
mat_path  = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code/epsilon_daily.mat';
out_root  = '/Users/turjo/Desktop/run_obs_eps_master';

% two time steps for robustness
dt_coarse = 1.0;   % days
dt_fine   = 0.5;   % days

% velocity-law parameters for export
PARAMS.D0_um = 1000;   % reference diameter (1 mm)
PARAMS.w0_md = 50;     % sinking speed at D0 (m d^-1)
PARAMS.beta  = 1.0;    % w ~ (D/D0)^beta
PARAMS.rho_p = 1.2;    % particle density (g cm^-3)

% ------------------------------------------------------------------------
% 1) LOAD ε(t) (robust)
% ------------------------------------------------------------------------
Sraw = load(mat_path);
[S, t_days, eps_series, span_str, rng_str] = grab_time_eps(Sraw);

% output dir
outdir = fullfile(out_root, ts_now());
if ~exist(outdir,'dir'), mkdir(outdir); end

fprintf('MAT file: %s\n', mat_path);
fprintf('Time span: %s\n', span_str);
fprintf('ε range  : %s\n', rng_str);
fprintf('Out dir : %s\n', outdir);

% ------------------------------------------------------------------------
% 2) RESAMPLE ε(t) to coarse/fine
% ------------------------------------------------------------------------
[tq_c, eps_c] = resample_eps(t_days, eps_series, dt_coarse);
[tq_f, eps_f] = resample_eps(t_days, eps_series, dt_fine);

% ------------------------------------------------------------------------
% 3) RUN MODEL (coarse & fine) — tolerate missing ε-hook
% ------------------------------------------------------------------------
[sim_c, used_hook_c] = run_sim_with_eps(tq_c, eps_c, dt_coarse);
[sim_f, used_hook_f] = run_sim_with_eps(tq_f, eps_f, dt_fine);

% ------------------------------------------------------------------------
% 4) EXTRACT PSD (absolute) and diameters
% ------------------------------------------------------------------------
[D_um, N_c] = get_psd_from_sim(sim_c, numel(tq_c));
assert(~isempty(D_um) && ~isempty(N_c), 'Could not locate PSD/D in outputs (coarse).');
[~,    N_f] = get_psd_from_sim(sim_f, numel(tq_f));
assert(~isempty(N_f), 'Could not locate PSD in outputs (fine).');

% ------------------------------------------------------------------------
% 5) EXPORT from PSD × w(D)
% ------------------------------------------------------------------------
mass_flux_c = compute_flux_from_psd(D_um, N_c, PARAMS);   % g cm^-2 d^-1
mass_flux_f = compute_flux_from_psd(D_um, N_f, PARAMS);   % g cm^-2 d^-1

% ------------------------------------------------------------------------
% 6) FIGURES
% ------------------------------------------------------------------------
make_overview_figure(tq_c, eps_c, D_um, N_c, mass_flux_c, used_hook_c, outdir);
make_mass_vs_settling_QA(tq_c, mass_flux_c, outdir);
make_dt_robust_plot(tq_c, mass_flux_c, tq_f, mass_flux_f, outdir);

fprintf('✓ All done. Outputs in: %s\n', outdir);

end % === main ===


% ========================================================================
% Helpers
% ========================================================================

function [S, t_days, eps_series, span_str, rng_str] = grab_time_eps(Sraw)
% Robustly find time and epsilon fields in (possibly nested) struct

% unwrap single top-level field if present
if isstruct(Sraw) && numel(fieldnames(Sraw))==1
    S = Sraw.(fieldnames(Sraw){1});
else
    S = Sraw;
end

% candidates (case-insensitive)
time_keys = {'mtime','time','t_days','days','t','datetime','datenum'};
eps_keys  = {'eps','epsilon','dissipation','epsilon_wkg','eps_series'};

[t_any, hit_t]   = lookup_field_recursive(S, time_keys);
[eps_any, hit_e] = lookup_field_recursive(S, eps_keys);

assert(~isempty(t_any),  'grab_time_eps: No time field found (looked for mtime/time/...).');
assert(~isempty(eps_any),'grab_time_eps: No epsilon field found (looked for eps/epsilon/...).');

% time to days-since-start (column)
t_any = t_any(:);
if any(t_any > 700000) || any(t_any > 1e5)     % datenum-like
    t_days = t_any - t_any(1);
else
    t_days = t_any - t_any(1);
end

% eps shape handling — collapse depth if present
E = squeeze(eps_any);
if size(E,1) == numel(t_days)
    eps_series = E(:,1);      % already aligned
elseif size(E,2) == numel(t_days)
    eps_series = E(1,:).';
else
    % If it's [nz x nt], average over depth
    if ndims(E)==2
        if size(E,2) == numel(t_days)
            eps_series = squeeze(nanmean(E,1)).';
        else
            eps_series = squeeze(nanmean(E,2));
        end
    else
        error('grab_time_eps: Unexpected epsilon shape.');
    end
end

span_str = sprintf('%.2f to %.2f days (N=%d)', t_days(1), t_days(end), numel(t_days));
rng_str  = sprintf('[%.2e, %.2e]', min(eps_series), max(eps_series));

fprintf('grab_time_eps: using time key "%s" (N=%d)\n', hit_t, numel(t_days));
fprintf('grab_time_eps: using epsilon key "%s"\n', hit_e);
end


function [val, hitkey] = lookup_field_recursive(S, keys)
% recursive case-insensitive search for any of KEYS inside struct S
val = []; hitkey = '';
if ~isstruct(S), return; end
fns = fieldnames(S);
% direct hits first
for i = 1:numel(keys)
    k = keys{i};
    j = find(strcmpi(fns, k), 1);
    if ~isempty(j)
        val = S.(fns{j});
        hitkey = fns{j};
        return;
    end
end
% recurse
for j = 1:numel(fns)
    v = S.(fns{j});
    if isstruct(v)
        [val, hitkey] = lookup_field_recursive(v, keys);
        if ~isempty(val), return; end
    end
end
end


function [tq, eps_q] = resample_eps(t_days, eps_series, dt_days)
tq = (t_days(1):dt_days:t_days(end)).';
eps_q = fillmissing(eps_series(:),'linear','EndValues','nearest');
eps_q = interp1(t_days, eps_q, tq, 'linear','extrap');
eps_q = smoothdata(eps_q,'movmedian',3);
end


function [sim, used_hook] = run_sim_with_eps(tq, eps_q, dt_days)
% Try to run CoagulationSimulation and wire epsilon if possible
sim = CoagulationSimulation();

% set dt if exposed
if isprop(sim,'config') && isprop(sim.config,'time_step_days')
    sim.config.time_step_days = dt_days;
elseif isprop(sim,'config') && isprop(sim.config,'dt_days')
    sim.config.dt_days = dt_days;
end

sim.initializeComponents();

used_hook = false;
% attempt common hooks
if ismethod(sim,'setEpsilon')
    try
        sim.setEpsilon(tq, eps_q);
        used_hook = true;
    catch
    end
elseif isprop(sim,'config')
    try
        sim.config.forcing.epsilon = [tq eps_q];
        used_hook = true;
    catch
    end
end

fprintf('ε forcing hook used: %d\n', used_hook);

% run (signature may be run() or run(t))
try
    sim.run(tq);
catch
    sim.run();
end
end


function [D_um, N_TxK] = get_psd_from_sim(sim, T_expected)
% Find diameters and PSD in common locations. Return [T x K].
D_um = [];
N_TxK = [];

% diameters
if isprop(sim,'grid')
    if isprop(sim.grid,'D_um')
        D_um = sim.grid.D_um(:);
    elseif isprop(sim.grid,'D')
        D_um = sim.grid.D(:) * 1e4; % cm -> µm
    end
end

% PSD
cands = {};
if isprop(sim,'outputs'), cands{end+1} = sim.outputs; end
if isprop(sim,'state'),   cands{end+1} = sim.state;   end
cands{end+1} = sim;

N = [];
for ii = 1:numel(cands)
    S = cands{ii};
    if isstruct(S)
        if isfield(S,'N'),  N = S.N;  break; end
        if isfield(S,'Nt'), N = S.Nt; break; end
    end
end
if isempty(N), return; end

N = squeeze(N);
% Make [T x K]
if size(N,1) == T_expected
    N_TxK = N;
elseif size(N,2) == T_expected
    N_TxK = N.';  % transpose
else
    % try to guess using min(abs(size - T_expected))
    [~,ix] = min([abs(size(N,1)-T_expected), abs(size(N,2)-T_expected)]);
    if ix==1, N_TxK = N;
    else,     N_TxK = N.'; end
end

% If diameters missing, synthesize
if isempty(D_um)
    K = size(N_TxK,2);
    D_um = logspace(log10(10), log10(1e4), K).';
end
end


function mass_flux = compute_flux_from_psd(D_um, N_TxK, P)
% Volume flux = sum_k N * volume_k * w_k  (cm^3 cm^-2 d^-1)
% Mass flux   = rho_p * Volume flux       (g cm^-2 d^-1)
D_um_row = D_um(:).';
D_cm     = D_um_row * 1e-4;                  % µm → cm
vol_cm3  = (pi/6) * (D_cm.^3);               % sphere
w_md     = P.w0_md * (D_um_row / P.D0_um).^P.beta;  % m d^-1
w_cm_d   = w_md * 100;                       % cm d^-1

mass_flux = P.rho_p * (N_TxK .* vol_cm3) * w_cm_d.'; % [T x 1]
mass_flux = make_col(mass_flux);
end


function make_overview_figure(t, eps_t, D_um, N_TxK, mass_flux, used_hook, outdir)
fig = figure('Color','w','Position',[100 100 1100 700]);

% ε(t)
subplot(2,2,1);
plot(t, eps_t, 'b-','LineWidth',1.1); grid on; set(gca,'YScale','log');
yline(1e-6,'k--','eref'); ylabel('\epsilon  (W kg^{-1})','Interpreter','tex');
xlabel('Time (days)'); title('\epsilon(t)','Interpreter','tex');

% Export time series
subplot(2,2,2);
plot(t, mass_flux, 'b-','LineWidth',1.2); grid on;
ylabel('Mass flux (g cm^{-2} d^{-1})','Interpreter','tex'); xlabel('Time (days)');
ttl = 'Export (from PSD \times w(D))';
if ~used_hook, ttl = [ttl '  —  ε hook NOT used']; end
title(ttl,'Interpreter','tex');

% PSD start vs end (absolute)
subplot(2,2,3);
loglog(D_um, max(N_TxK(1,:),realmin), 'b-','LineWidth',1.1); hold on;
loglog(D_um, max(N_TxK(end,:),realmin),'r-','LineWidth',1.1); grid on;
xlabel('Diameter  (\mum)','Interpreter','tex'); ylabel('#   cm^{-3}','Interpreter','tex');
legend(sprintf('t=%.1f d',t(1)), sprintf('t=%.1f d',t(end)),'Location','southwest');
title('PSD at start vs end (absolute)','Interpreter','tex');

% Heatmap log10 N
subplot(2,2,4);
imagesc(t, log10(D_um), log10(max(N_TxK,realmin)));
axis xy; cb = colorbar; ylabel(cb,'log_{10}  #  cm^{-3}','Interpreter','tex');
xlabel('Time (days)'); ylabel('log_{10}  D  (\mum)','Interpreter','tex');
title('log_{10}  PSD  (time \times size)','Interpreter','tex');

save_fig(fig, outdir, 'observed_eps.png');
end


function make_mass_vs_settling_QA(t, mass_flux, outdir)
% Here "settling integral" = cumulative integral of mass flux
fig = figure('Color','w','Position',[100 100 1200 380]);
mass_int = cumtrapz(t, max(mass_flux,0));
plot(t, mass_int, 'r--','LineWidth',1.6); grid on; hold on;
% keep a thin blue zero line for visual consistency
plot(t, zeros(size(t)), 'b-','LineWidth',1.1);
ylabel('Mass / Loss (model units)','Interpreter','tex');
xlabel('Time (days)'); title('Mass vs Settling (QA)','Interpreter','tex');
legend(' \int  settling  dt', 'Mass lost (proxy)','Location','northwest');
save_fig(fig, outdir, 'mass_vs_settling.png');
end


function make_dt_robust_plot(t1, F1, t2, F2, outdir)
% Interpolate fine onto coarse for a max-diff %
F2_on_1 = interp1(t2, F2, t1, 'linear','extrap');
den     = max(max(abs([F1(:);F2_on_1(:)])), eps);
mdiff   = 100*max(abs(F1(:)-F2_on_1(:)))/den;

fig = figure('Color','w','Position',[100 100 1200 380]);
plot(t1, F1, 'b-','LineWidth',1.2); hold on;
plot(t2, F2, 'r--','LineWidth',1.2); grid on;
xlabel('Time (days)'); ylabel('Mass flux (g cm^{-2} d^{-1})');
title(sprintf('\\Delta t robustness (max diff = %.2f%%)', mdiff),'Interpreter','tex');
legend('\Delta t = 1 d','\Delta t = 0.5 d','Location','best');
save_fig(fig, outdir, 'dt_robustness.png');
end


function save_fig(fig, outdir, name)
fn = fullfile(outdir, name);
try
    saveas(fig, fn);
catch
    warning('Could not save figure: %s', fn);
end
end


function y = make_col(x)
y = x(:);
end


function s = ts_now()
c = clock;
s = sprintf('%04d%02d%02d_%02d%02d%02.0f', c(1),c(2),c(3),c(4),c(5),floor(c(6)));
end