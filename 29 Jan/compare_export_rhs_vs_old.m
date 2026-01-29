clear; close all; clc;

% --- EDIT THIS ---
phase1_dir = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/26 Jan/results_phase1_20260127_152841';
assert(isfolder(phase1_dir), 'Not a folder: %s', phase1_dir);

coderoot = fileparts(phase1_dir);
addpath(genpath(coderoot));

% Load saved run
S1 = load(fullfile(phase1_dir,'cfg_snapshot.mat')); cfg = S1.cfg;
S2 = load(fullfile(phase1_dir,'out_baseline.mat')); out = S2.out;

t = out.time(:);
Y = out.concentrations;

% Rebuild sim (no ODE solve)
sim = CoagulationSimulation(cfg);
sim.buildRHS();

% --------- biovolume weights (same space as your M(t)) ----------
gridS = sim.rhs.getSectionGrid();
try
    w = gridS.getBinVolumes();  % cm^3 per particle per bin
catch
    if ismethod(gridS,'getConservedRadii')
        r_cm = gridS.getConservedRadii();
    else
        r_cm = gridS.getFractalRadii();
    end
    w = (4/3)*pi*(r_cm(:).^3);
end
w = w(:);

Ns = numel(w);
Nz = cfg.getNumLayers();
dz_cm = cfg.dz * 100;

intCol = @(vflat) sum(sum(reshape(vflat,[Ns Nz]).*repmat(w,1,Nz))) * dz_cm;

% --------- compute M(t), dMdt_fd, and RHS-integrated terms ----------
nt = numel(t);
M         = zeros(nt,1);
dMdt_sum  = zeros(nt,1);   % full sum of terms
PP_rhs    = zeros(nt,1);   % integrated dv_pp

for it = 1:nt
    v = Y(it,:)';
    M(it) = intCol(v);

    T = sim.rhs.decomposeTerms(t(it), v);

    % safe fill
    if ~isfield(T,'dv_lin')   || isempty(T.dv_lin),    T.dv_lin    = zeros(size(v)); end
    if ~isfield(T,'dv_coag')  || isempty(T.dv_coag),   T.dv_coag   = zeros(size(v)); end
    if ~isfield(T,'dv_disagg')|| isempty(T.dv_disagg), T.dv_disagg = zeros(size(v)); end
    if ~isfield(T,'dv_other') || isempty(T.dv_other),  T.dv_other  = zeros(size(v)); end
    if ~isfield(T,'dv_pp')    || isempty(T.dv_pp),     T.dv_pp     = zeros(size(v)); end

    dMdt_sum(it) = intCol(T.dv_lin + T.dv_coag + T.dv_disagg + T.dv_other + T.dv_pp);
    PP_rhs(it)   = intCol(T.dv_pp);
end

dMdt_fd = nan(nt,1);
dMdt_fd(2:end) = diff(M)./diff(t);

% --------- KEY: RHS-consistent "export-like sink" ----------
Export_rhs = PP_rhs - dMdt_sum;

% RHS-consistent simple residual (should be tiny like full closure)
RES_simple_rhs = dMdt_fd - (PP_rhs - Export_rhs);  % equals dMdt_fd - dMdt_sum

ok = isfinite(dMdt_fd);
fprintf('median |RES_simple_rhs| = %.3e\n', median(abs(RES_simple_rhs(ok))));
fprintf('max    |RES_simple_rhs| = %.3e\n', max(abs(RES_simple_rhs(ok))));

% --------- OLD export (your diagnostic) ----------
Export_old = nan(nt,1);
src = 'none';

if isfield(out,'output_data') && ~isempty(out.output_data)
    od = out.output_data;

    % try common names used in your code base
    if isfield(od,'total_flux_best_cm3cm2d') && ~isempty(od.total_flux_best_cm3cm2d)
        Export_old = od.total_flux_best_cm3cm2d(:);
        src = 'out.output_data.total_flux_best_cm3cm2d';
    elseif isfield(od,'bottom_flux_cm3cm2d') && ~isempty(od.bottom_flux_cm3cm2d)
        Export_old = od.bottom_flux_cm3cm2d(:);
        src = 'out.output_data.bottom_flux_cm3cm2d';
    end
end

fprintf('Export_old source: %s\n', src);

% --------- PLOTS ----------
outpng = fullfile(phase1_dir,'compare_export_rhs_vs_old.png');

fig = figure('Color','w','Position',[80 80 1200 700]);
tiledlayout(2,1,'Padding','compact','TileSpacing','compact');

nexttile; hold on; grid on;
plot(t, Export_rhs, 'LineWidth', 2);
if any(isfinite(Export_old))
    plot(t, Export_old, '--', 'LineWidth', 2);
    legend({'Export\_rhs (PP\_rhs - dMdt\_sum)','Export\_old (diagnostic)'},'Location','best');
else
    legend({'Export\_rhs (PP\_rhs - dMdt\_sum)'},'Location','best');
end
xlabel('t [d]'); ylabel('rate (same units as M(t) per day)');
title('Export comparison in inventory space');

nexttile; hold on; grid on;
if any(isfinite(Export_old))
    plot(t, Export_old - Export_rhs, 'LineWidth', 2);
    yline(0,'k-');
    title('Difference: Export\_old - Export\_rhs');
    xlabel('t [d]'); ylabel('difference');
else
    text(0.1,0.5,'Export\_old not found in out.output\_data','Units','normalized');
    axis off;
end

exportgraphics(fig, outpng, 'Resolution', 200);
close(fig);

fprintf('Saved: %s\n', outpng);