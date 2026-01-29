% ============================================================
% run_phase2_tendencies.m  (UPDATED: grouped figures)
%
% Phase-2 (detective) plots using a SAVED Phase-1 folder
% - Loads cfg + out from the Phase-1 results folder
% - Rebuilds sim (no need to re-run ODE)
% - Produces dQ/dt vs size plots (total + by process) at 3 depths
% - Times: 19–25 days
%
% OUTPUTS (saved in a new folder inside phase1 folder):
%   phase2_tendencies_YYYYMMDD_HHMMSS/
%     dQdt_tot_all_depths.png
%     dQdt_coag_all_depths.png
%     dQdt_lin_all_depths.png
%     dQdt_disagg_all_depths.png
%     dQdt_other_all_depths.png
%     dQdt_pp_all_depths.png
%
% IMPORTANT:
% - For growth_mode='pp', the "PP-as-growth" term is usually inside dv_other.
%   So: look at dv_other for PP forcing by size (in your current setup).
% ============================================================

clear; close all; clc;

phase1_dir = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/26 Jan/results_phase1_20260127_173801';
assert(isfolder(phase1_dir), 'Not a folder: %s', phase1_dir);

% Add code root (parent folder of results_phase1_*)
coderoot = fileparts(phase1_dir);
addpath(genpath(coderoot));

% Load cfg + out
S1 = load(fullfile(phase1_dir,'cfg_snapshot.mat')); cfg = S1.cfg;
S2 = load(fullfile(phase1_dir,'out_baseline.mat')); out = S2.out;

t = out.time(:);
Y = out.concentrations;

% Rebuild simulation objects (no ODE solve)
sim = CoagulationSimulation(cfg);

% Build betas + RHS (exists in your 26 Jan version)
sim.buildRHS();

% Output folder for Phase-2 plots (inside Phase-1 folder)
stamp  = datestr(now,'yyyymmdd_HHMMSS');
outdir = fullfile(phase1_dir, ['phase2_tendencies_' stamp]);
if ~exist(outdir,'dir'), mkdir(outdir); end

% ------------------------------------------------------------
% Choose depths and times
% ------------------------------------------------------------
z = cfg.getZ(); z = z(:);

% requested depths (m)
zreq_list = [10, 50, 100];
if max(z) < 100
    zreq_list(3) = max(z); % deepest if column is shallow
end

kz_list = zeros(size(zreq_list));
zpick_list = zeros(size(zreq_list));
for i = 1:numel(zreq_list)
    [~,kz_list(i)] = min(abs(z - zreq_list(i)));
    zpick_list(i)  = z(kz_list(i));
end

% times (days)
time_list = 19:25;

% ------------------------------------------------------------
% Size axis (D in um)
% ------------------------------------------------------------
D_um = [];
try
    if isfield(out,'output_data') && isfield(out.output_data,'diam_i') && ~isempty(out.output_data.diam_i)
        D_um = out.output_data.diam_i(:) * 1e4; % cm -> um
    end
catch
end

if isempty(D_um)
    % fallback using radii
    gridS = sim.rhs.getSectionGrid();
    if ismethod(gridS,'getConservedRadii')
        r_cm = gridS.getConservedRadii();
    else
        r_cm = gridS.getFractalRadii();
    end
    D_um = 2*r_cm(:)*1e4;
end

[D_sort, isort] = sort(D_um(:));

% Ns/Nz reshape params
gridS = sim.rhs.getSectionGrid();
if ismethod(gridS,'getConservedRadii')
    r_cm = gridS.getConservedRadii(); r_cm = r_cm(:);
else
    r_cm = gridS.getFractalRadii();   r_cm = r_cm(:);
end
Ns = numel(r_cm);
Nz = cfg.getNumLayers();

% ------------------------------------------------------------
% Grouped tendency figures (ONE figure per term, 3 panels)
% ------------------------------------------------------------
make_term_plot_grouped(sim, t, Y, Ns, Nz, kz_list, zpick_list, D_sort, isort, time_list, ...
    'dv_tot',  'total dQ/dt',  outdir);

make_term_plot_grouped(sim, t, Y, Ns, Nz, kz_list, zpick_list, D_sort, isort, time_list, ...
    'dv_coag', 'coagulation tendency', outdir);

make_term_plot_grouped(sim, t, Y, Ns, Nz, kz_list, zpick_list, D_sort, isort, time_list, ...
    'dv_lin',  'linear/transport+sinking tendency', outdir);

make_term_plot_grouped(sim, t, Y, Ns, Nz, kz_list, zpick_list, D_sort, isort, time_list, ...
    'dv_disagg','disaggregation tendency', outdir);

make_term_plot_grouped(sim, t, Y, Ns, Nz, kz_list, zpick_list, D_sort, isort, time_list, ...
    'dv_other','other tendency (often includes PP-as-growth)', outdir);

make_term_plot_grouped(sim, t, Y, Ns, Nz, kz_list, zpick_list, D_sort, isort, time_list, ...
    'dv_pp',   'PP source tendency (explicit PP only)', outdir);

fprintf('\nPHASE-2 grouped tendency plots saved in:\n  %s\n', outdir);

% ============================================================
% Functions
% ============================================================
function make_term_plot_grouped(sim, t, Y, Ns, Nz, kz_list, zpick_list, D_sort, isort, time_list, field, pretty, outdir)

fig = figure('Color','w','Position',[100 100 1500 500]);
set(fig,'Toolbar','none'); set(fig,'Menubar','none');

tl = tiledlayout(1, numel(kz_list), 'TileSpacing','compact', 'Padding','compact');

for iz = 1:numel(kz_list)
    kz    = kz_list(iz);
    zpick = zpick_list(iz);

    nexttile; hold on;

    for jt = 1:numel(time_list)
        treq = time_list(jt);
        [~,itp] = min(abs(t - treq));

        vflat = Y(itp,:)';
        T = sim.rhs.decomposeTerms(t(itp), vflat);

        if ~isfield(T,field) || isempty(T.(field))
            continue;
        end

        dv  = T.(field);
        dv2 = reshape(dv, [Ns Nz]);
        y   = dv2(:,kz);
        y   = y(:);
        y   = y(isort);

        plot(D_sort, y, 'LineWidth', 1.4);
    end

    yline(0,'k-');
    grid on;
    set(gca,'XScale','log');

    xlabel('D [\mum]');
    ylabel('tendency (state units per day)');
    title(sprintf('z = %.1f m', zpick));

    xlim([max(1, min(D_sort(D_sort>0))) max(D_sort)]);
end

title(tl, sprintf('%s (lines = days %d–%d)', pretty, time_list(1), time_list(end)));
legend(string(time_list), 'Location','eastoutside');

fname = sprintf('dQdt_%s_all_depths.png', sanitize_name(field));
exportgraphics(fig, fullfile(outdir, fname), 'Resolution', 200);
close(fig);

end

function s = sanitize_name(field)
% make filename-friendly
s = regexprep(field,'^dv_','');          % dv_tot -> tot
s = regexprep(s,'[^a-zA-Z0-9]+','_');    % safe
end