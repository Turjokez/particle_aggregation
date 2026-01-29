function [outdir, dbg] = run_Diagnostics_FINAL(cfg)
% run_Diagnostics_FINAL(cfg)
% Diagnostics suite (single entry point).
%
% UPDATE (BIOVOLUME BUDGET FIX):
% - DOES NOT require cfg.diag_budget_space (your SimulationConfig does not have it).
% - Forces MASS BALANCE diagnostics to run in BIOVOLUME space (Vbin weights).
% - Uses OutputGenerator bottom flux if available.
% - Correctly handles cm3/cm2/day vs cm3/m2/day flux units.

clc;
dbg = struct();

% ==========================================================
% 0) EPSILON FORCING (OPTIONAL, SAFE)
% ==========================================================
use_eps_file = true;
eps_file = 'epsilon_daily.mat';

has_series = isprop(cfg,'epsilon_time') && isprop(cfg,'epsilon_series') && ...
            ~isempty(cfg.epsilon_time) && ~isempty(cfg.epsilon_series);
has_fun    = isprop(cfg,'eps_fun') && ~isempty(cfg.eps_fun);

if use_eps_file && ~has_series && ~has_fun
    try
        A = load(eps_file);
        if isfield(A,'S'), S = A.S; else, S = A; end

        assert(isfield(S,'mtime') && isfield(S,'z') && isfield(S,'eps'), ...
            'epsilon_daily.mat must contain S.mtime, S.z, S.eps');

        t_days  = S.mtime(:) - S.mtime(1);

        % --- get model z safely ---
        z_model = [];
        try
            z_model = cfg.getZ();
        catch
        end
        if isempty(z_model)
            if isprop(cfg,'use_column') && cfg.use_column
                Nz_tmp = cfg.getNumLayers();
                z_model = ((0:Nz_tmp-1) + 0.5) * cfg.dz;
            else
                z_model = 0;
            end
        end
        z_model = z_model(:);

        z_data  = S.z(:);

        E  = S.eps;
        Nt = numel(t_days);

        % ensure E is [Nz_data x Nt]
        if size(E,2) ~= Nt && size(E,1) == Nt
            E = E.';
        end
        assert(size(E,2) == Nt, 'eps must have Nt columns after orientation fix');

        % interpolate in z to model z
        E_tz    = E.'; % Nt x Nz_data
        E_model = interp1(z_data, E_tz.', z_model, 'linear', 'extrap').'; % Nt x Nz_model

        cfg.epsilon_time   = t_days;
        cfg.epsilon_series = E_model;

        if isprop(cfg,'epsilon_const'),      cfg.epsilon_const = []; end
        if isprop(cfg,'epsilon_interface'),  cfg.epsilon_interface = 'series'; end

        fprintf('Loaded epsilon forcing from %s\n', eps_file);
        fprintf('  time: [%d]  eps(model): [%s]\n', numel(cfg.epsilon_time), mat2str(size(cfg.epsilon_series)));

        dbg.eps_loaded = true;
        dbg.eps_time_len = numel(cfg.epsilon_time);
        dbg.eps_series_size = size(cfg.epsilon_series);

    catch ME
        fprintf('WARNING: epsilon file load failed (%s). Continuing without it.\n', ME.message);
        dbg.eps_loaded = false;
        dbg.eps_error = ME.message;
    end
else
    if has_series && isprop(cfg,'epsilon_interface'), cfg.epsilon_interface = 'series'; end
    if has_fun    && isprop(cfg,'epsilon_interface'), cfg.epsilon_interface = 'fun';    end
end

% ==========================================================
% 1) RUN SIM
% ==========================================================
sim = CoagulationSimulation(cfg);
out = sim.run();

t = out.time(:);
Y = out.concentrations;

% ==========================================================
% 1.5) ENSURE out.output_data EXISTS (SAFE)
% ==========================================================
try
    has_od = isfield(out,'output_data') && ~isempty(out.output_data);
catch
    has_od = false;
end

if ~has_od
    try
        out.output_data = OutputGenerator.spectraAndFluxes(t, Y, sim.grid, cfg);
        dbg.output_data_created = true;
    catch ME
        dbg.output_data_created = false;
        dbg.output_data_error = ME.message;
        fprintf('WARNING: could not compute output_data (%s). Continuing.\n', ME.message);
    end
else
    dbg.output_data_created = false;
end

% ==========================================================
% 2) OUTPUT FOLDER
% ==========================================================
stamp  = datestr(now,'yyyymmdd_HHMMSS');
outdir = fullfile(pwd, ['results_diag_' stamp]);
if ~exist(outdir,'dir'), mkdir(outdir); end

save(fullfile(outdir,'cfg_snapshot.mat'),'cfg');
save(fullfile(outdir,'out_baseline.mat'),'out','-v7.3');

fprintf('\nSaved outputs to:\n  %s\n\n', outdir);

figdir_spectra = fullfile(outdir,'fig_spectra');
figdir_tend    = fullfile(outdir,'fig_tendencies');
if ~exist(figdir_spectra,'dir'), mkdir(figdir_spectra); end
if ~exist(figdir_tend,'dir'),    mkdir(figdir_tend);    end

% ==========================================================
% 3) GRID + GEOMETRY + WEIGHT
% ==========================================================
gridS = sim.rhs.getSectionGrid();

if ismethod(gridS,'getConservedRadii')
    r_cm = gridS.getConservedRadii(); r_cm = r_cm(:);
else
    r_cm = gridS.getFractalRadii();   r_cm = r_cm(:);
end

Ns = [];
try
    if isprop(cfg,'n_sections') && ~isempty(cfg.n_sections)
        Ns = cfg.n_sections;
    end
catch
end
if isempty(Ns), Ns = numel(r_cm); end
if isempty(Ns)
    error('Could not determine Ns. Add cfg.n_sections or make sure grid has radii.');
end

is_col = isprop(cfg,'use_column') && cfg.use_column;
if is_col
    Nz = cfg.getNumLayers();
else
    Nz = 1;
end

dz_cm = cfg.dz * 100;

% weights
w_num = ones(Ns,1);

w_vbin = [];
try
    if ismethod(gridS,'getBinVolumes')
        w_vbin = gridS.getBinVolumes(); w_vbin = w_vbin(:);
    end
catch
end
if isempty(w_vbin)
    try
        w_vbin = Disaggregation.getBinVolumes_cm3(gridS); w_vbin = w_vbin(:);
    catch
    end
end
if isempty(w_vbin)
    w_vbin = (4/3)*pi*(r_cm(:).^3);
end

% auto weight (for coag conservation checks only)
[w_auto, wname_auto, wdbg] = getConservedWeight_AUTO(cfg, sim, gridS, r_cm, Ns, Nz, dz_cm, t, Y);
dbg.weight_debug = wdbg;

% ==========================================================
% IMPORTANT: FORCE biovolume-space diagnostics here
% (no cfg.diag_budget_space required)
% ==========================================================
w_budget = w_vbin;
wname_budget = "Vbin(budget)";

dbg.is_col = is_col;
dbg.Nz = Nz;
dbg.Ns = Ns;
dbg.weight_auto_name = wname_auto;
dbg.weight_budget_name = wname_budget;
dbg.weight_budget_min = min(w_budget);
dbg.weight_budget_max = max(w_budget);

% ==========================================================
% 3.5) DIAMETERS (for thresholds)
% ==========================================================
D_um = [];
try
    if isfield(out,'output_data') && isfield(out.output_data,'diam_i') && ~isempty(out.output_data.diam_i)
        D_um = out.output_data.diam_i(:) * 1e4; % cm -> um
        dbg.D_um_source = 'out.output_data.diam_i';
    end
catch
end
if isempty(D_um)
    D_um = 2 * r_cm * 1e4;
    dbg.D_um_source = 'r_cm (fallback)';
end
dbg.D_um_min = min(D_um);
dbg.D_um_max = max(D_um);

% ==========================================================
% 4) EXPORT FLUX (prefer OutputGenerator)
% ==========================================================
sink_rate_vec = [];
if isfield(sim.operators,'sink_rate') && ~isempty(sim.operators.sink_rate)
    sink_rate_vec = sim.operators.sink_rate(:);
end

fprintf('DEBUG: is_col=%d Nz=%d | sink_rate_vec=%s | weight_auto=%s | budget=%s\n', ...
    is_col, Nz, mat2str(size(sink_rate_vec)), wname_auto, wname_budget);
dbg.sink_rate_vec_size = size(sink_rate_vec);

Fsect_export = [];
Fsect_name   = '';

try
    if isfield(out,'output_data') && ~isempty(out.output_data)
        od = out.output_data;

        if isfield(od,'bottom_fluxsect_cm3cm2d') && ~isempty(od.bottom_fluxsect_cm3cm2d)
            Fsect_export = od.bottom_fluxsect_cm3cm2d; % [nt x Ns]
            Fsect_name   = 'bottom_fluxsect_cm3cm2d';
        elseif isfield(od,'bottom_fluxsect_cm3m2d') && ~isempty(od.bottom_fluxsect_cm3m2d)
            Fsect_export = od.bottom_fluxsect_cm3m2d;  % [nt x Ns]
            Fsect_name   = 'bottom_fluxsect_cm3m2d';
        end
    end
catch
end

dbg.export_flux_source = Fsect_name;
dbg.export_sinkrate_used_only_if_no_output_flux = isempty(Fsect_export);

% ==========================================================
% 5) SIZE CLASS THRESHOLDS
% ==========================================================
thr_small = 500;
thr_mid_A = 2000;
thr_mid_B = 1000;

idxA_small = find(D_um < thr_small);
idxA_med   = find(D_um >= thr_small & D_um < thr_mid_A);
idxA_large = find(D_um >= thr_mid_A);

idxB_small = find(D_um < thr_small);
idxB_med   = find(D_um >= thr_small & D_um < thr_mid_B);
idxB_large = find(D_um >= thr_mid_B);

dbg.nA = [numel(idxA_small) numel(idxA_med) numel(idxA_large)];
dbg.nB = [numel(idxB_small) numel(idxB_med) numel(idxB_large)];

% ==========================================================
% 6) MASS BALANCE (BIOVOLUME SPACE)
% ==========================================================
inventory_mode = "integrated";
nt = numel(t);

M   = zeros(nt,1);
EX  = zeros(nt,1);

dMdt_lin    = zeros(nt,1);
dMdt_coag   = zeros(nt,1);
dMdt_disagg = zeros(nt,1);
dMdt_other  = zeros(nt,1);
dMdt_pp     = zeros(nt,1);

intTerm = @(dv) integrateMassProxy(dv, w_budget, Ns, Nz, dz_cm, inventory_mode);

use_pp = false;
try
    use_pp = isprop(cfg,'enable_pp') && cfg.enable_pp;
catch
end

for it = 1:nt
    vflat = Y(it,:)';

    % inventory in BV space
    M(it) = intTerm(vflat);

    T = sim.rhs.decomposeTerms(t(it), vflat);

    if isfield(T,'dv_lin')    && ~isempty(T.dv_lin),    dMdt_lin(it)    = intTerm(T.dv_lin);    end
    if isfield(T,'dv_coag')   && ~isempty(T.dv_coag),   dMdt_coag(it)   = intTerm(T.dv_coag);   end
    if isfield(T,'dv_disagg') && ~isempty(T.dv_disagg), dMdt_disagg(it) = intTerm(T.dv_disagg); end
    if isfield(T,'dv_other')  && ~isempty(T.dv_other),  dMdt_other(it)  = intTerm(T.dv_other);  end

    if use_pp && isfield(T,'dv_pp') && ~isempty(T.dv_pp)
        dMdt_pp(it) = intTerm(T.dv_pp);
    else
        dMdt_pp(it) = 0;
    end

    % export EX in the SAME BV space
    if ~isempty(Fsect_export)
        F = max(Fsect_export(it,:), 0);

        if strcmp(Fsect_name, 'bottom_fluxsect_cm3m2d')
            EX(it) = sum(F) * 1e-4;   % m^2 -> cm^2
        else
            EX(it) = sum(F);          % already cm^3/cm^2/day
        end
    else
        % fallback (only if no OutputGenerator flux exists)
        if is_col && Nz > 1 && ~isempty(sink_rate_vec)
            N2 = reshape(vflat,[Ns Nz]);
            vbot = N2(:,end);
            ex_num = vbot .* max(sink_rate_vec,0);
            EX(it) = sum(ex_num .* w_budget(:));
        else
            EX(it) = 0;
        end
    end
end

dMdt_fd = nan(nt,1);
dMdt_fd(2:end) = diff(M)./diff(t);

RES_simple = dMdt_fd - (dMdt_pp - EX);
RES_full   = dMdt_fd - (dMdt_lin + dMdt_coag + dMdt_disagg + dMdt_other + dMdt_pp);

res_ok  = isfinite(RES_simple) & ((1:nt)' > 1);
fprintf('Budget residual (simple, %s): median=%.3e max=%.3e\n', ...
    wname_budget, median(abs(RES_simple(res_ok))), max(abs(RES_simple(res_ok))));

res_ok2  = isfinite(RES_full) & ((1:nt)' > 1);
fprintf('Budget residual (FULL terms, %s): median=%.3e max=%.3e\n', ...
    wname_budget, median(abs(RES_full(res_ok2))), max(abs(RES_full(res_ok2))));

fprintf('DEBUG: median(|dMdt_disagg|)=%.3e | max(|dMdt_disagg|)=%.3e\n', ...
    median(abs(dMdt_disagg(isfinite(dMdt_disagg)))), max(abs(dMdt_disagg(isfinite(dMdt_disagg)))));

dbg.M = M; dbg.EX = EX;
dbg.dMdt_fd = dMdt_fd;
dbg.dMdt_lin = dMdt_lin;
dbg.dMdt_coag = dMdt_coag;
dbg.dMdt_disagg = dMdt_disagg;
dbg.dMdt_other = dMdt_other;
dbg.dMdt_pp = dMdt_pp;
dbg.RES_simple = RES_simple;
dbg.RES_full = RES_full;

fig = figure('Color','w'); turnOffToolbars(fig);
plot(t,RES_simple,'LineWidth',1.5); grid on;
xlabel('t [d]'); ylabel('Residual');
title(sprintf('Budget residual (simple, %s): dM/dt - (PP - EX)', wname_budget));
safeExportPNG(fig, fullfile(outdir,'mass_balance_residual_simple.png'));
if ~isempty(fig) && isvalid(fig), close(fig); end

fig = figure('Color','w'); turnOffToolbars(fig);
plot(t,RES_full,'LineWidth',1.5); grid on;
xlabel('t [d]'); ylabel('Residual');
title(sprintf('Budget residual (FULL, %s): dM/dt - \\int(sum terms)', wname_budget));
safeExportPNG(fig, fullfile(outdir,'mass_balance_residual_full.png'));
if ~isempty(fig) && isvalid(fig), close(fig); end

% ==========================================================
% 7) EXPORT SIZE FRACTIONS (from flux if available)
% ==========================================================
if ~isempty(Fsect_export)
    [FsB,FmB,FlB,tot_export] = exportFractions_FromFlux(Fsect_export, idxB_small, idxB_med, idxB_large);
    [FsA,FmA,FlA,~]          = exportFractions_FromFlux(Fsect_export, idxA_small, idxA_med, idxA_large);
else
    [FsB,FmB,FlB,tot_export] = exportFractions(Y, Ns, Nz, w_auto, sink_rate_vec, idxB_small, idxB_med, idxB_large);
    [FsA,FmA,FlA,~]          = exportFractions(Y, Ns, Nz, w_auto, sink_rate_vec, idxA_small, idxA_med, idxA_large);
end

dbg.Fs_500_1000 = FsB; dbg.Fm_500_1000 = FmB; dbg.Fl_500_1000 = FlB;
dbg.Fs_500_2000 = FsA; dbg.Fm_500_2000 = FmA; dbg.Fl_500_2000 = FlA;
dbg.tot_export = tot_export;

fig = figure('Position',[80 80 1200 650], 'Color','w'); turnOffToolbars(fig);
plot(t,FsB,'LineWidth',2.0); hold on;
plot(t,FmB,'LineWidth',2.0);
plot(t,FlB,'LineWidth',2.0);
grid on; ylim([0 1]);
xlabel('t [d]'); ylabel('fraction');
title('Export size fractions (500/1000)');
legend({'<500','500-1000','>=1000'},'Location','best');
safeExportPNG(fig, fullfile(outdir,'export_sizefractions_500_1000.png'));
if ~isempty(fig) && isvalid(fig), close(fig); end

fig = figure('Position',[80 80 1200 650], 'Color','w'); turnOffToolbars(fig);
plot(t,FsA,'LineWidth',2.0); hold on;
plot(t,FmA,'LineWidth',2.0);
plot(t,FlA,'LineWidth',2.0);
grid on; ylim([0 1]);
xlabel('t [d]'); ylabel('fraction');
title('Export size fractions (500/2000)');
legend({'<500','500-2000','>=2000'},'Location','best');
safeExportPNG(fig, fullfile(outdir,'export_sizefractions_500_2000.png'));
if ~isempty(fig) && isvalid(fig), close(fig); end

% ==========================================================
% 7.5) SIZE SPECTRA (ONLY if column)
% ==========================================================
do_spectra = true;

if do_spectra && is_col && isfield(out,'output_data') && isfield(out.output_data,'Y3')
    od = out.output_data;

    depth_list_m = [10 50 100];
    time_list_d  = [0 5 10 19 20 21 22 23 24 25];

    zc  = od.z(:);
    Dum = od.diam_i(:) * 1e4;

    [Dum_sort, isort] = sort(Dum);

    fig = figure('Color','w','Position',[80 80 1200 350]); turnOffToolbars(fig);
    tl = tiledlayout(1, numel(depth_list_m), 'Padding','compact', 'TileSpacing','compact');

    for iz = 1:numel(depth_list_m)
        zreq = depth_list_m(iz);
        [~,kz] = min(abs(zc - zreq));
        zpick = zc(kz);

        nexttile; hold on;
        for jt = 1:numel(time_list_d)
            treq = time_list_d(jt);
            [~,itp] = min(abs(t - treq));

            Nlayer = squeeze(od.Y3(itp, kz, :));
            Nplot  = max(Nlayer(:), 0);
            Nplot  = Nplot(isort);

            loglog(Dum_sort, Nplot, 'LineWidth', 1.0);
        end

        grid on;
        xlabel('D [\mum]'); ylabel('N');
        title(sprintf('z=%.1f m', zpick));
        set(gca,'XLim',[max(1,min(Dum_sort)) max(Dum_sort)]);
    end

    title(tl, 'Size spectra (lines = selected times)');
    legend(string(time_list_d), 'Location','eastoutside');
    safeExportPNG(fig, fullfile(figdir_spectra, 'sizespectra_all_depths.png'));
    if ~isempty(fig) && isvalid(fig), close(fig); end
end

% ==========================================================
% 7.6) TENDENCIES vs SIZE (ONLY if column)
% ==========================================================
do_tendencies = true;

if do_tendencies && is_col
    depth_list_m = [10 50 100];
    time_list_d  = [19 20 21 22 23 24 25];

    Dum = D_um(:);
    [Dum_sort, isort] = sort(Dum);

    zc = cfg.getZ(); zc = zc(:);

    term_names = {'dv_tot','dv_pp','dv_coag','dv_disagg','dv_lin','dv_other'};
    pretty     = {'total','PP','coag','disagg','linear_transport_sink','other'};

    for tn = 1:numel(term_names)
        field = term_names{tn};

        fig = figure('Color','w','Position',[80 80 1200 350]); turnOffToolbars(fig);
        tl = tiledlayout(1, numel(depth_list_m), 'Padding','compact', 'TileSpacing','compact');

        for iz = 1:numel(depth_list_m)
            zreq = depth_list_m(iz);
            [~,kz] = min(abs(zc - zreq));
            zpick = zc(kz);

            nexttile; hold on;

            for jt = 1:numel(time_list_d)
                treq = time_list_d(jt);
                [~,itp] = min(abs(t - treq));

                vflat = Y(itp,:)';
                T = sim.rhs.decomposeTerms(t(itp), vflat);

                if ~isfield(T,field) || isempty(T.(field))
                    continue;
                end

                dv  = T.(field);
                dv2 = reshape(dv, [Ns Nz]);
                dvk = dv2(:, kz);

                y = dvk(:);
                y = y(isort);
                plot(Dum_sort, y, 'LineWidth', 1.0);
            end

            grid on;
            set(gca,'XScale','log');
            xlabel('D [\mum]'); ylabel('tendency');
            title(sprintf('z=%.1f m', zpick));
        end

        title(tl, sprintf('dQ/dt vs size | term=%s (lines = times)', pretty{tn}));
        legend(string(time_list_d), 'Location','eastoutside');

        fname = sprintf('dQdt_%s_all_depths.png', pretty{tn});
        safeExportPNG(fig, fullfile(figdir_tend, fname));
        if ~isempty(fig) && isvalid(fig), close(fig); end
    end
end

% ==========================================================
% 8) CLOSURE CHECK
% ==========================================================
[cl_rel, cl_abs] = closureCheck(sim, t, Y);
fprintf('Closure check: max rel inf = %.3e | max abs inf = %.3e\n', cl_rel, cl_abs);
dbg.closure_max_relinf = cl_rel;
dbg.closure_max_absinf = cl_abs;

% ==========================================================
% 9) SAVE
% ==========================================================
save(fullfile(outdir,'dbg.mat'),'dbg');
fprintf('\nDONE. Folder:\n  %s\n', outdir);

end

% ================= helpers =================

function [w, wname, dbg] = getConservedWeight_AUTO(cfg, sim, gridS, r_cm, Ns, Nz, dz_cm, t, Y)
dbg = struct();

wOnes = ones(Ns,1);

wVbin = [];
try
    if ismethod(gridS,'getBinVolumes')
        wVbin = gridS.getBinVolumes();
        wVbin = wVbin(:);
    end
catch
end
if isempty(wVbin)
    try
        wVbin = Disaggregation.getBinVolumes_cm3(gridS);
        wVbin = wVbin(:);
    catch
    end
end
if isempty(wVbin)
    wVbin = (4/3)*pi*(r_cm(:).^3);
end

wAv = [];
try
    if isprop(sim.grid,'av_vol') && ~isempty(sim.grid.av_vol)
        wAv = sim.grid.av_vol(:);
    end
catch
end
if isempty(wAv)
    wAv = wVbin;
end

t0 = t(1);
v0 = Y(1,:).';

T  = sim.rhs.decomposeTerms(t0, v0);
dv = T.dv_coag;

I = @(x,wbin) integrateMassProxy(x, wbin, Ns, Nz, dz_cm, "integrated");

T0_ones = I(v0, wOnes);
dT_ones = I(dv, wOnes);

T0_vbin = I(v0, wVbin);
dT_vbin = I(dv, wVbin);

T0_av   = I(v0, wAv);
dT_av   = I(dv, wAv);

leakOnes = abs(dT_ones) / max(abs(T0_ones), 1e-30);
leakVbin = abs(dT_vbin) / max(abs(T0_vbin), 1e-30);
leakAv   = abs(dT_av)   / max(abs(T0_av),   1e-30);

dbg.leakOnes = leakOnes;
dbg.leakVbin = leakVbin;
dbg.leakAv   = leakAv;

[~,ix] = min([leakOnes, leakVbin, leakAv]);

if ix==1
    w = wOnes; wname = "ones(auto)";
elseif ix==2
    w = wVbin; wname = "Vbin(auto)";
else
    w = wAv;   wname = "av_vol(auto)";
end

fprintf('DEBUG(weight-auto): leak ones=%.3e | leak Vbin=%.3e | leak av=%.3e -> using %s\n', ...
    leakOnes, leakVbin, leakAv, wname);
end

function val = integrateMassProxy(dv, w, Ns, Nz, dz_cm, inventory_mode)
if isempty(dv), val = 0; return; end
if Nz==1
    val = sum(dv(:) .* w(:));
else
    N2 = reshape(dv(:), [Ns Nz]);
    W2 = repmat(w(:), 1, Nz);
    if inventory_mode == "surface"
        val = sum(N2(:,1) .* w(:));
    else
        val = sum(sum(N2 .* W2)) * dz_cm;
    end
end
end

function [Fs,Fm,Fl,tot_export] = exportFractions_FromFlux(Fsect, idxS, idxM, idxL)
nt = size(Fsect,1);
Fs = nan(nt,1); Fm = nan(nt,1); Fl = nan(nt,1);
tot_export = nan(nt,1);

for it = 1:nt
    F = max(Fsect(it,:),0);
    tot = sum(F);
    tot_export(it) = tot;

    if tot > 0 && isfinite(tot)
        Fs(it) = sum(F(idxS))/tot;
        Fm(it) = sum(F(idxM))/tot;
        Fl(it) = sum(F(idxL))/tot;
    end
end
end

function [Fs,Fm,Fl,tot_export] = exportFractions(Y, Ns, Nz, w, sink_rate_vec, idxS, idxM, idxL)
nt = size(Y,1);
Fs = nan(nt,1); Fm = nan(nt,1); Fl = nan(nt,1);
tot_export = nan(nt,1);

if Nz==1 || isempty(sink_rate_vec)
    return;
end

for it = 1:nt
    vflat = Y(it,:)';
    N2 = reshape(vflat,[Ns Nz]);
    vbot = N2(:,end);

    ex_num  = vbot .* max(sink_rate_vec,0);
    ex_w    = ex_num .* w(:);

    tot = sum(ex_w);
    tot_export(it) = tot;

    if tot > 0 && isfinite(tot)
        Fs(it) = sum(ex_w(idxS))/tot;
        Fm(it) = sum(ex_w(idxM))/tot;
        Fl(it) = sum(ex_w(idxL))/tot;
    end
end
end

function [max_rel, max_abs] = closureCheck(sim, t, Y)
max_rel = 0; max_abs = 0;
ii = unique(round(linspace(1, numel(t), 8)));

for k = 1:numel(ii)
    it = ii(k);
    v = Y(it,:).';
    terms = sim.rhs.decomposeTerms(t(it), v);

    if ~isfield(terms,'dv_other') || isempty(terms.dv_other)
        terms.dv_other = zeros(size(v));
    end
    if ~isfield(terms,'dv_pp') || isempty(terms.dv_pp)
        terms.dv_pp = zeros(size(v));
    end
    if ~isfield(terms,'dv_disagg') || isempty(terms.dv_disagg)
        terms.dv_disagg = zeros(size(v));
    end

    dv_sum = terms.dv_lin + terms.dv_coag + terms.dv_pp + terms.dv_disagg + terms.dv_other;
    dv_err = terms.dv_tot - dv_sum;

    abs_err = norm(dv_err, inf);
    scale   = max(1e-30, norm(terms.dv_tot, inf));
    rel_err = abs_err / scale;

    max_abs = max(max_abs, abs_err);
    max_rel = max(max_rel, rel_err);
end
end

function turnOffToolbars(fig)
try, set(fig,'Toolbar','none'); end %#ok<TRYNC>
try, set(fig,'Menubar','none'); end %#ok<TRYNC>
end

function safeExportPNG(fig, filename)
if isempty(fig) || ~isvalid(fig)
    fprintf('WARNING: figure invalid before export: %s\n', filename);
    return;
end

try
    drawnow;

    try
        exportgraphics(fig, filename, 'Resolution', 200);
        return;
    catch
    end

    set(fig,'PaperPositionMode','auto');
    print(fig, filename, '-dpng', '-r200');

catch ME
    fprintf('WARNING: failed to export %s (%s)\n', filename, ME.message);
end
end