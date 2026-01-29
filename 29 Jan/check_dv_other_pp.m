clear; close all; clc;

% --- EDIT THIS ---
phase1_dir = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/26 Jan/Result/Baseline';
assert(isfolder(phase1_dir), 'Not a folder: %s', phase1_dir);

coderoot = fileparts(phase1_dir);
addpath(genpath(coderoot));

S1 = load(fullfile(phase1_dir,'cfg_snapshot.mat')); cfg = S1.cfg;
S2 = load(fullfile(phase1_dir,'out_baseline.mat')); out = S2.out;

t = out.time(:);
Y = out.concentrations;

sim = CoagulationSimulation(cfg);
sim.buildRHS();

% --- weights for biovolume budget (same idea as your diagnostics) ---
gridS = sim.rhs.getSectionGrid();
try
    w = gridS.getBinVolumes();   % best if exists
catch
    % fallback (sphere)
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

% --- loop times and integrate tendencies ---
nt = numel(t);
dMdt_tot   = zeros(nt,1);
dMdt_sum   = zeros(nt,1);
dMdt_other = zeros(nt,1);
dMdt_pp    = zeros(nt,1);

for it = 1:nt
    v = Y(it,:)';
    T = sim.rhs.decomposeTerms(t(it), v);

    % safe fill
    if ~isfield(T,'dv_other') || isempty(T.dv_other), T.dv_other = zeros(size(v)); end
    if ~isfield(T,'dv_pp')    || isempty(T.dv_pp),    T.dv_pp    = zeros(size(v)); end
    if ~isfield(T,'dv_lin')   || isempty(T.dv_lin),   T.dv_lin   = zeros(size(v)); end
    if ~isfield(T,'dv_coag')  || isempty(T.dv_coag),  T.dv_coag  = zeros(size(v)); end
    if ~isfield(T,'dv_disagg')|| isempty(T.dv_disagg),T.dv_disagg= zeros(size(v)); end

    dMdt_tot(it)   = intCol(T.dv_tot);
    dMdt_other(it) = intCol(T.dv_other);
    dMdt_pp(it)    = intCol(T.dv_pp);

    dMdt_sum(it)   = intCol(T.dv_lin + T.dv_coag + T.dv_disagg + T.dv_other + T.dv_pp);
end

% --- print summary numbers ---
ii = (2:nt)'; % skip first point
fprintf('\n=== Integrated tendency summary (biovolume space) ===\n');
fprintf('median |dMdt_tot|   = %.3e\n', median(abs(dMdt_tot(ii))));
fprintf('median |dMdt_sum|   = %.3e\n', median(abs(dMdt_sum(ii))));
fprintf('median |dMdt_other| = %.3e\n', median(abs(dMdt_other(ii))));
fprintf('median |dMdt_pp|    = %.3e\n', median(abs(dMdt_pp(ii))));

fprintf('\nmedian |(tot - sum)| = %.3e\n', median(abs(dMdt_tot(ii)-dMdt_sum(ii))));
fprintf('max    |(tot - sum)| = %.3e\n', max(abs(dMdt_tot(ii)-dMdt_sum(ii))));

% --- very direct PP test ---
if median(abs(dMdt_pp(ii))) < 1e-20 && median(abs(dMdt_other(ii))) > 0
    fprintf('\nNOTE: dv_pp ~ 0 but dv_other is nonzero -> PP/growth is likely inside dv_other.\n');
end