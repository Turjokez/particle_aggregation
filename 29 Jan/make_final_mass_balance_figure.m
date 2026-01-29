% ============================================================
% make_final_mass_balance_figure.m
% FINAL mass-balance figure from a SAVED Phase-1 folder
% - simple residual
% - full residual (includes dv_other etc.)
% - cumulative curves
%
% Output:
%   <phase1_dir>/diagnostics_mass_balance_FINAL.png
% ============================================================

clear; close all; clc;

% --- UPDATE THIS ---
phase1_dir = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/26 Jan/results_phase1_20260127_190210';
assert(isfolder(phase1_dir), 'Not a folder: %s', phase1_dir);

% add code root so classes are on path
coderoot = fileparts(phase1_dir);
addpath(genpath(coderoot));

% output
out_png = fullfile(phase1_dir, 'diagnostics_mass_balance_FINAL.png');

% ------------------------------------------------------------
% 1) Try to load dbg.mat (fast path)
% ------------------------------------------------------------
dbg = [];
dbg_file = fullfile(phase1_dir, 'dbg.mat');
if isfile(dbg_file)
    try
        A = load(dbg_file);
        if isfield(A,'dbg')
            dbg = A.dbg;
        end
    catch
        warning('dbg.mat exists but could not be read. Recomputing diagnostics.');
        dbg = [];
    end
end

% Always load cfg + out (needed for time axis)
S1 = load(fullfile(phase1_dir,'cfg_snapshot.mat')); cfg = S1.cfg;
S2 = load(fullfile(phase1_dir,'out_baseline.mat')); out = S2.out;

t  = out.time(:);
nt = numel(t);

% ------------------------------------------------------------
% 2) If dbg missing OR if PP/EX not present, recompute
% ------------------------------------------------------------
need_compute = isempty(dbg) || ~isfield(dbg,'M') || ~isfield(dbg,'EX') || ~isfield(dbg,'dMdt_fd') || ...
               ~isfield(dbg,'dMdt_pp') || ~isfield(dbg,'RES_simple') || ~isfield(dbg,'RES_full');

if need_compute
    fprintf('dbg.mat missing/incomplete -> computing mass balance terms...\n');

    Y = out.concentrations;
    sim = CoagulationSimulation(cfg);
    if ismethod(sim,'buildRHS'), sim.buildRHS(); end

    % ---- grid/weights for BIOVOLUME inventory ----
    gridS = sim.rhs.getSectionGrid();

    if ismethod(gridS,'getConservedRadii')
        r_cm = gridS.getConservedRadii(); r_cm = r_cm(:);
    else
        r_cm = gridS.getFractalRadii();   r_cm = r_cm(:);
    end
    Ns = numel(r_cm);
    Nz = cfg.getNumLayers();
    dz_cm = cfg.dz * 100;

    % biovolume weight per bin (cm^3)
    w_vbin = [];
    try
        if ismethod(gridS,'getBinVolumes')
            w_vbin = gridS.getBinVolumes(); w_vbin = w_vbin(:);
        end
    catch
    end
    if isempty(w_vbin)
        w_vbin = (4/3)*pi*(r_cm.^3);
    end

    % bottom export flux by section if available in output_data
    Fsect_export = [];
    F_is_m2 = false;
    try
        od = out.output_data;
        if isfield(od,'bottom_fluxsect_cm3cm2d') && ~isempty(od.bottom_fluxsect_cm3cm2d)
            Fsect_export = od.bottom_fluxsect_cm3cm2d; % [nt x Ns]
            F_is_m2 = false;
        elseif isfield(od,'bottom_fluxsect_cm3m2d') && ~isempty(od.bottom_fluxsect_cm3m2d)
            Fsect_export = od.bottom_fluxsect_cm3m2d;  % [nt x Ns]
            F_is_m2 = true;
        end
    catch
    end

    % helper integration (column-integrated BV)
    intTerm = @(vflat) integrateMassProxy(vflat, w_vbin, Ns, Nz, dz_cm);

    M           = zeros(nt,1);
    EX          = zeros(nt,1);

    dMdt_fd     = nan(nt,1);
    dMdt_lin    = zeros(nt,1);
    dMdt_coag   = zeros(nt,1);
    dMdt_disagg = zeros(nt,1);
    dMdt_other  = zeros(nt,1);
    dMdt_pp     = zeros(nt,1);

    for it = 1:nt
        vflat = Y(it,:)';
        M(it) = intTerm(vflat);

        T = sim.rhs.decomposeTerms(t(it), vflat);

        if isfield(T,'dv_lin')    && ~isempty(T.dv_lin),    dMdt_lin(it)    = intTerm(T.dv_lin);    end
        if isfield(T,'dv_coag')   && ~isempty(T.dv_coag),   dMdt_coag(it)   = intTerm(T.dv_coag);   end
        if isfield(T,'dv_disagg') && ~isempty(T.dv_disagg), dMdt_disagg(it) = intTerm(T.dv_disagg); end
        if isfield(T,'dv_other')  && ~isempty(T.dv_other),  dMdt_other(it)  = intTerm(T.dv_other);  end

        % --------------------------
        % PP term (this is the fix)
        % --------------------------
        if isfield(T,'dv_pp') && ~isempty(T.dv_pp)
            dMdt_pp(it) = intTerm(T.dv_pp);
        else
            dMdt_pp(it) = 0;
        end

        % export in same units: cm^3/cm^2/day
        if ~isempty(Fsect_export)
            F = max(Fsect_export(it,:),0);
            if F_is_m2
                EX(it) = sum(F) * 1e-4; % m^2 -> cm^2
            else
                EX(it) = sum(F);
            end
        end
    end

    if nt >= 2
        dMdt_fd(2:end) = diff(M)./diff(t);
        dMdt_fd(1) = dMdt_fd(2);
    else
        dMdt_fd(:) = 0;
    end

    RES_simple = dMdt_fd - (dMdt_pp - EX);
    RES_full   = dMdt_fd - (dMdt_lin + dMdt_coag + dMdt_disagg + dMdt_other + dMdt_pp);

    dbg = struct();
    dbg.M = M; dbg.EX = EX;
    dbg.dMdt_fd = dMdt_fd;
    dbg.dMdt_lin = dMdt_lin;
    dbg.dMdt_coag = dMdt_coag;
    dbg.dMdt_disagg = dMdt_disagg;
    dbg.dMdt_other = dMdt_other;
    dbg.dMdt_pp = dMdt_pp;
    dbg.RES_simple = RES_simple;
    dbg.RES_full = RES_full;

    try
        save(fullfile(phase1_dir,'dbg_massbalance_recompute.mat'),'dbg');
    catch
    end
end

% ------------------------------------------------------------
% 3) Build cumulative curves
% ------------------------------------------------------------
M = dbg.M(:);
EX = dbg.EX(:);

dMdt_fd    = dbg.dMdt_fd(:);
dMdt_pp    = dbg.dMdt_pp(:);
dMdt_lin   = getfield_or(dbg,'dMdt_lin', zeros(nt,1));
dMdt_coag  = getfield_or(dbg,'dMdt_coag', zeros(nt,1));
dMdt_disagg= getfield_or(dbg,'dMdt_disagg', zeros(nt,1));
dMdt_other = getfield_or(dbg,'dMdt_other', zeros(nt,1));

net_simple = (dMdt_pp - EX);
net_full   = (dMdt_lin + dMdt_coag + dMdt_disagg + dMdt_other + dMdt_pp);

Cum = @(x) cumtrapz(t, x);
dM = M - M(1);

I_simple = Cum(net_simple);
I_full   = Cum(net_full);

RES_simple = dbg.RES_simple(:);
RES_full   = dbg.RES_full(:);

% ------------------------------------------------------------
% 4) Plot one FINAL figure (multi-panel)
% ------------------------------------------------------------
fig = figure('Color','w','Position',[60 60 1200 900]);
set(fig,'Toolbar','none'); set(fig,'Menubar','none');

tl = tiledlayout(3,2,'Padding','compact','TileSpacing','compact');

% (1) Inventory
nexttile;
plot(t, M, 'LineWidth', 2); grid on;
xlabel('t [d]'); ylabel('M(t) [cm^3/cm^2]');
title('Column inventory (biovolume)');

% (2) Rates
nexttile;
hold on;
plot(t, dMdt_fd,  'LineWidth', 2);
plot(t, dMdt_pp,  'LineWidth', 2);
plot(t, EX,       'LineWidth', 2);
plot(t, net_full, 'LineWidth', 1.5);
grid on;
xlabel('t [d]'); ylabel('Rate [cm^3/cm^2/d]');
title('Rates: dM/dt, PP, Export, Net(full)');
legend({'dM/dt (FD)','PP (dv\_pp)','Export','Net(full)'}, 'Location','best');

% (3) Simple residual
nexttile;
plot(t, RES_simple, 'LineWidth', 2); grid on;
xlabel('t [d]'); ylabel('Residual');
title('Residual (simple): dM/dt - (PP - Export)');
yline(0,'k-');

% (4) Full residual
nexttile;
plot(t, RES_full, 'LineWidth', 2); grid on;
xlabel('t [d]'); ylabel('Residual');
title('Residual (full): dM/dt - (sum of terms)');
yline(0,'k-');

% (5) Cumulative
nexttile([1 2]);
hold on;
plot(t, dM,       'LineWidth', 2);
plot(t, I_simple, 'LineWidth', 2);
plot(t, I_full,   'LineWidth', 2);
grid on;
xlabel('t [d]'); ylabel('Cumulative [cm^3/cm^2]');
title('Cumulative: \DeltaM, \int(PP-EX)dt, \int(Net(full))dt');
legend({'\DeltaM','\int(PP-EX)dt','\int(Net(full))dt'}, 'Location','best');

title(tl, sprintf('FINAL Mass Balance Summary | %s', phase1_dir), 'Interpreter','none');

exportgraphics(fig, out_png, 'Resolution', 200);
close(fig);

fprintf('Saved: %s\n', out_png);

% ============================================================
% helpers
% ============================================================
function val = integrateMassProxy(vflat, w, Ns, Nz, dz_cm)
if Nz == 1
    val = sum(vflat(:) .* w(:));
else
    N2 = reshape(vflat(:), [Ns Nz]);
    W2 = repmat(w(:), 1, Nz);
    val = sum(sum(N2 .* W2)) * dz_cm;
end
end

function x = getfield_or(S, name, fallback)
if isfield(S,name) && ~isempty(S.(name))
    x = S.(name);
else
    x = fallback;
end
end