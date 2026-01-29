function S = diag_biovolume_budget_run(sim, out)
%DIAG_BIOVOLUME_BUDGET_RUN  Simple biovolume budget check:
%   M(t)  vs  M0 + ∫(PP - EX)dt
%
% Works with your current 26 Jan output_data fields (preferred):
%   od.state_units
%   od.total_mass_best_cm3cm2              (BEST inventory, cm^3/cm^2)
%   od.total_flux_best_cm3cm2d             (BEST export, cm^3/cm^2/day)
%   od.bottom_total_flux_weighted_cm3cm2d
%   od.bottom_fluxsect_weighted_cm3cm2d
%   od.diam_i
%
% Also supports older fallback fields (kept):
%   od.column_total_inventory_cm3cm2
%   od.bottom_total_flux_cm3cm2d
%   od.total_mass, od.total_flux, od.fluxsect
%
% Usage:
%   S = diag_biovolume_budget_run(csim, out);

t   = out.time(:);
Y   = out.concentrations;
od  = out.output_data;
cfg = sim.config;

Ns  = cfg.n_sections;
Nz  = cfg.getNumLayers();
dzcm = cfg.dz * 100;

% -------------------------------
% 0) state_units
% -------------------------------
state_units = "unknown";
if isfield(od,'state_units') && ~isempty(od.state_units)
    state_units = string(od.state_units);
end
fprintf('state_units = %s\n', state_units);

% -------------------------------
% 1) Inventory M(t) in cm3/cm2
% -------------------------------
% PREFERRED: use "best" cm^2 inventory when available
M_bv = [];

if isfield(od,'total_mass_best_cm3cm2') && ~isempty(od.total_mass_best_cm3cm2)
    M_bv = od.total_mass_best_cm3cm2(:);
    fprintf('Using od.total_mass_best_cm3cm2 for inventory (mean=%g)\n', mean(M_bv,'omitnan'));
elseif isfield(od,'column_total_inventory_weighted_cm3cm2') && ~isempty(od.column_total_inventory_weighted_cm3cm2)
    M_bv = od.column_total_inventory_weighted_cm3cm2(:);
    fprintf('Using od.column_total_inventory_weighted_cm3cm2 for inventory (mean=%g)\n', mean(M_bv,'omitnan'));
elseif isfield(od,'column_total_inventory_cm3cm2') && ~isempty(od.column_total_inventory_cm3cm2)
    M_bv = od.column_total_inventory_cm3cm2(:);
    fprintf('Using od.column_total_inventory_cm3cm2 for inventory (mean=%g)\n', mean(M_bv,'omitnan'));
end

% --- OLD LOGIC (kept, but now used only to build a cross-check) ---
% If state is cm3/cm3 (biovolume concentration), then:
%   column inventory per area = sum_over_bins,sum_over_z( Y * dzcm )
if state_units == "cm3"
    if Nz == 1
        M_from_state_A = sum(Y, 2) * dzcm;                 % assume Y is cm3/cm3
        M_from_state_B = sum(Y, 2);                        % assume Y is already cm3/cm2 (no dz)
    else
        Nt = numel(t);
        M_from_state_A = zeros(Nt,1);
        M_from_state_B = zeros(Nt,1);
        for it = 1:Nt
            A = reshape(Y(it,:).', [Ns Nz]);
            M_from_state_A(it) = sum(A(:)) * dzcm;         % assume cm3/cm3
            M_from_state_B(it) = sum(A(:));                % assume already integrated
        end
    end
else
    M_from_state_A = [];
    M_from_state_B = [];
end

% If M_bv missing, fall back to legacy od.total_mass (kept)
if isempty(M_bv)
    if isfield(od,'total_mass') && ~isempty(od.total_mass)
        M_bv = od.total_mass(:);
        warning('No cm^2 inventory field found. Using od.total_mass as inventory (check units!).');
    else
        error('Cannot build inventory: no od.total_mass_best_cm3cm2 / column_total_inventory_* and no od.total_mass found.');
    end
end

% Choose which "from state" check is more consistent (A uses dzcm, B does not)
M_from_state = [];
if ~isempty(M_from_state_A) && ~isempty(M_from_state_B)
    errA = max(abs(M_bv - M_from_state_A));
    errB = max(abs(M_bv - M_from_state_B));
    if errA <= errB
        M_from_state = M_from_state_A;
        fprintf('Inventory check (state * dzcm): max|M_bv - M_from_state| = %.3e\n', errA);
    else
        M_from_state = M_from_state_B;
        fprintf('Inventory check (state no dz):  max|M_bv - M_from_state| = %.3e\n', errB);
    end
elseif ~isempty(M_from_state_A)
    fprintf('Inventory check: max|M_bv - M_from_state| = %.3e\n', max(abs(M_bv - M_from_state_A)));
    M_from_state = M_from_state_A;
end

% -------------------------------
% Inventory used in budget
% -------------------------------
% IMPORTANT: do NOT "force" to M_from_state by default anymore.
% Use the output_data BEST fields as truth.
M_budget = M_bv;

% (optional) still report mismatch
if ~isempty(M_from_state)
    inv_mismatch = max(abs(M_bv - M_from_state));
    fprintf('Inventory mismatch report: max|M_bv - M_from_state| = %.3e\n', inv_mismatch);
end

% -------------------------------
% 2) Export EX(t) in cm3/cm2/day
% -------------------------------
EX_bv = [];

% PREFERRED: use "best" cm^2 export when available
if isfield(od,'total_flux_best_cm3cm2d') && ~isempty(od.total_flux_best_cm3cm2d)
    EX_bv = od.total_flux_best_cm3cm2d(:);
    fprintf('Using od.total_flux_best_cm3cm2d for EX (mean=%g)\n', mean(EX_bv,'omitnan'));
elseif isfield(od,'bottom_total_flux_weighted_cm3cm2d') && ~isempty(od.bottom_total_flux_weighted_cm3cm2d)
    EX_bv = od.bottom_total_flux_weighted_cm3cm2d(:);
    fprintf('Using od.bottom_total_flux_weighted_cm3cm2d for EX (mean=%g)\n', mean(EX_bv,'omitnan'));
elseif isfield(od,'bottom_total_flux_cm3cm2d') && ~isempty(od.bottom_total_flux_cm3cm2d)
    EX_bv = od.bottom_total_flux_cm3cm2d(:);
    fprintf('Using od.bottom_total_flux_cm3cm2d for EX (mean=%g)\n', mean(EX_bv,'omitnan'));
end

% Cross-check from per-bin export (if present)
if isfield(od,'bottom_fluxsect_weighted_cm3cm2d') && ~isempty(od.bottom_fluxsect_weighted_cm3cm2d)
    EX_from_bins = sum(od.bottom_fluxsect_weighted_cm3cm2d, 2);
    fprintf('Check: mean(sum(bottom_fluxsect_weighted_cm3cm2d)) = %g\n', mean(EX_from_bins,'omitnan'));
    if isempty(EX_bv)
        EX_bv = EX_from_bins;
        fprintf('Using sum(od.bottom_fluxsect_weighted_cm3cm2d) for EX.\n');
    else
        fprintf('Check: max|EX_bv - EX_from_bins| = %g\n', max(abs(EX_bv - EX_from_bins)));
    end
end

% --- OLD LOGIC (kept): od.total_flux / od.fluxsect ---
if isempty(EX_bv)
    if isfield(od,'total_flux') && ~isempty(od.total_flux)
        EX_bv = od.total_flux(:);
        fprintf('Using legacy od.total_flux for EX (mean=%g)\n', mean(EX_bv,'omitnan'));
    else
        if ~isfield(od,'fluxsect') || isempty(od.fluxsect)
            error('No export flux found (no total_flux_best_cm3cm2d, no bottom_total_flux*, no bottom_fluxsect_weighted_cm3cm2d, no od.total_flux, no od.fluxsect).');
        end

        if state_units == "cm3"
            EX_bv = sum(od.fluxsect, 2);
            warning('Using legacy od.fluxsect as already-biovolume (state_units=="cm3"). Check this is correct.');
        else
            Vbin = sim.grid.av_vol(:).';
            EX_bv = sum(od.fluxsect .* Vbin, 2);
            warning('state_units not "cm3": assuming od.fluxsect is number-flux and converting with Vbin.');
        end
    end
end

% -------------------------------
% 3) PP(t) in cm3/cm2/day from RHS dv_pp
% -------------------------------
Nt = numel(t);
PP_bv = zeros(Nt,1);

% preallocate
PP_from_dvpp_A = nan(Nt,1);
PP_from_dvpp_B = nan(Nt,1);

for it = 1:Nt
    v = Y(it,:).';
    T = sim.rhs.decomposeTerms(t(it), v);

    dv_pp = T.dv_pp(:);  % state/day

    if state_units == "cm3"
        % (A) cm3/cm3/day  -> multiply by dzcm
        % (B) cm3/cm2/day  -> do NOT multiply by dzcm
        PP_from_dvpp_A(it) = sum(dv_pp) * dzcm;
        PP_from_dvpp_B(it) = sum(dv_pp);
    else
        % dv_pp assumed #/cm3/day -> multiply by Vbin
        Vbin = sim.grid.av_vol(:);
        PP_bv(it) = sum(dv_pp .* Vbin) * dzcm;
    end
end

if state_units == "cm3"
    PP_bv = PP_from_dvpp_A(:); % temporary
end

% -------------------------------
% 3b) Implied PP from inventory tendency:  PP = dM/dt + EX
% -------------------------------
dMdt = gradient(M_budget, t);      % cm3/cm2/day
PP_implied = dMdt + EX_bv;         % cm3/cm2/day

if state_units == "cm3"
    errA = max(abs(PP_from_dvpp_A(:) - PP_implied(:)));
    errB = max(abs(PP_from_dvpp_B(:) - PP_implied(:)));
    if errB < errA
        PP_bv = PP_from_dvpp_B(:);
        fprintf('PP scaling choice: using dv_pp with NO dz (closer to implied).\n');
    else
        PP_bv = PP_from_dvpp_A(:);
        fprintf('PP scaling choice: using dv_pp * dzcm (closer to implied).\n');
    end
end

% -------------------------------
% 3c) PP used in budget
% -------------------------------
PP_budget = PP_bv; % default
if isprop(cfg,'enable_pp') && ~cfg.enable_pp
    fprintf('FORCING PP: cfg.enable_pp=false -> using PP_implied for simple budget.\n');
    PP_budget = PP_implied;
end

% -------------------------------
% 4) Budget compare
% -------------------------------
M_pred  = M_budget(1) + cumtrapz(t, PP_budget - EX_bv);
resid_M = M_budget - M_pred;

fprintf('max|resid| = %.3e\n', max(abs(resid_M)));
fprintf('rel max|resid| = %.3e\n', max(abs(resid_M))/max(abs(M_budget) + 1e-30));

dM_amp = max(M_budget) - min(M_budget);
fprintf('rel max|resid| (vs max inventory change) = %.3e\n', max(abs(resid_M))/max(dM_amp, 1e-30));

% -------------------------------
% 5) Plots
% -------------------------------
figure;
plot(t, M_budget, 'LineWidth',2); hold on; grid on;
plot(t, M_pred,'--','LineWidth',2);
legend('M (model)','M0 + \int(PP-EX)dt','Location','best');
xlabel('time (days)'); ylabel('cm^3/cm^2');
title('Inventory: model vs simple budget');

figure;
plot(t, resid_M, 'LineWidth',2); grid on;
xlabel('time (days)'); ylabel('cm^3/cm^2');
title('Simple budget residual (should be small)');

figure;
plot(t, PP_budget, 'LineWidth',2); hold on; grid on;
plot(t, EX_bv, 'LineWidth',2);
legend('PP (used)','EX','Location','best');
xlabel('time (days)'); ylabel('cm^3/cm^2/day');
title('Biovolume PP and export');

figure;
plot(t, PP_bv, 'LineWidth',2); hold on; grid on;
plot(t, PP_implied, '--', 'LineWidth',2);
legend('PP from dv\_pp','PP implied (dM/dt + EX)','Location','best');
xlabel('time (days)'); ylabel('cm^3/cm^2/day');
title('PP comparison (dv\_pp vs implied)');

% -------------------------------
% 6) Return struct
% -------------------------------
S = struct();
S.t = t;

S.M_bv = M_bv;
S.M_from_state = M_from_state;
S.M_budget = M_budget;

S.dMdt = dMdt;

S.PP_bv = PP_bv;
S.PP_budget = PP_budget;

S.EX_bv = EX_bv;
S.PP_implied = PP_implied;

S.M_pred = M_pred;
S.resid = resid_M;
S.state_units = state_units;

end