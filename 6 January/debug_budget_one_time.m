function debug_budget_one_time(simObj, t_query)
%DEBUG_BUDGET_ONE_TIME  One-time column budget diagnostic at a chosen time.
%
% Prints:
%   (A) Model budget using RHS decomposition:
%         dM/dt (RHS) = Int(dv_tot)
%                    = Int(dv_lin + dv_pp + dv_coag + dv_disagg)
%
%   (B) FD estimate from stored inventory M(t) (reference only)
%
% Notes:
% - Export is bottom export flux [cm^3 m^-2 d^-1], positive downward.
% - dv_* are tendencies in the SAME units used for inventory density
%   (biovolume concentration tendency), so column integration is unweighted.

    if nargin < 2
        error('Usage: debug_budget_one_time(simObj, t_query)');
    end
    if isempty(simObj) || ~isa(simObj,'CoagulationSimulation')
        error('debug_budget_one_time expects a CoagulationSimulation object.');
    end
    if ~ismethod(simObj.rhs,'decomposeTerms')
        error('simObj.rhs must implement decomposeTerms(t,v).');
    end

    out = simObj.result;
    cfg = simObj.config;
    grd = simObj.grid;

    if ~isfield(out,'time') || ~isfield(out,'concentrations')
        error('simObj.result must contain fields: time, concentrations.');
    end

    t = out.time(:);

    % Ensure output_data exists (needs M and export)
    if ~isfield(out,'output_data') || isempty(out.output_data)
        out.output_data = OutputGenerator.spectraAndFluxes(out.time, out.concentrations, grd, cfg);
        simObj.result.output_data = out.output_data; % optional: store back
    end
    od = out.output_data;

    % Closest saved time
    [~, it] = min(abs(t - t_query));
    tt = t(it);

    % Geometry
    Ns    = cfg.n_sections;
    Nz    = cfg.getNumLayers();
    dz_cm = cfg.dz * 100;

    % Inventory and export
    if ~isfield(od,'column_total_inventory_cm3m2') || isempty(od.column_total_inventory_cm3m2)
        error('output_data.column_total_inventory_cm3m2 missing.');
    end
    if ~isfield(od,'bottom_total_flux_cm3m2d') || isempty(od.bottom_total_flux_cm3m2d)
        error('output_data.bottom_total_flux_cm3m2d missing.');
    end

    M = od.column_total_inventory_cm3m2(:);
    F_export = od.bottom_total_flux_cm3m2d(:); % +down

    % FD (reference only)
    dMdt_fd_all = gradient(M, t);
    dMdt_fd = dMdt_fd_all(it);

    % State
    v = out.concentrations(it,:).';

    % Decompose RHS
    terms = simObj.rhs.decomposeTerms(tt, v);

    % Column integrals [cm^3 m^-2 d^-1]
    S_tot    = integrateColumnRate_local(terms.dv_tot,    Ns, Nz, dz_cm);
    S_lin    = integrateColumnRate_local(terms.dv_lin,    Ns, Nz, dz_cm);
    S_pp     = integrateColumnRate_local(terms.dv_pp,     Ns, Nz, dz_cm);
    S_coag   = integrateColumnRate_local(terms.dv_coag,   Ns, Nz, dz_cm);
    S_disagg = integrateColumnRate_local(terms.dv_disagg, Ns, Nz, dz_cm);

    S_sum = S_lin + S_pp + S_coag + S_disagg;

    % Print
    fprintf('\n=== Column budget @ t = %.4f d (requested %.2f) ===\n', tt, t_query);

    fprintf('Inventory M(t)                 = %+ .6e  [cm^3 m^-2]\n', M(it));
    fprintf('Export (bottom, +down)         = %+ .6e  [cm^3 m^-2 d^-1]\n', F_export(it));
    fprintf('\n');

    fprintf('dM/dt (RHS = Int(dv_tot))      = %+ .6e  [cm^3 m^-2 d^-1]\n', S_tot);
    fprintf('dM/dt (FD @ t=%.2f, ref only)  = %+ .6e  [cm^3 m^-2 d^-1]\n', tt, dMdt_fd);
    fprintf('FD mismatch: RHS - FD          = %+ .6e\n', S_tot - dMdt_fd);
    fprintf('\n');

    fprintf('Int(dv_tot)                    = %+ .6e\n', S_tot);
    fprintf('Int(dv_lin)                    = %+ .6e\n', S_lin);
    fprintf('Int(dv_pp)                     = %+ .6e\n', S_pp);
    fprintf('Int(dv_coag)                   = %+ .6e\n', S_coag);
    fprintf('Int(dv_disagg)                 = %+ .6e\n', S_disagg);
    fprintf('Sum(lin+pp+coag+disagg)         = %+ .6e\n', S_sum);
    fprintf('\n');

    fprintf('Closure (RHS): tot - sum(terms) = %+ .6e\n', S_tot - S_sum);

    % "Rough" sanity:
    % If dv_lin contains more than sinking (e.g., includes growth), this won't be ~0.
    fprintf('Sinking sanity (rough): Int(dv_lin) + Export = %+ .6e\n', S_lin + F_export(it));
    fprintf('===============================================\n');


    % --------------------------------------------------------------
    % Local helper (kept inside file to avoid dependency)
    % --------------------------------------------------------------
    function rate_cm3m2d = integrateColumnRate_local(dv, Ns, Nz, dz_cm)
        dv2 = reshape(dv(:), [Ns, Nz]);     % Ns x Nz
        layer_sum = sum(dv2, 1);            % 1 x Nz
        rate_cm3cm2d = sum(layer_sum) * dz_cm;
        rate_cm3m2d  = rate_cm3cm2d * 1e4;  % cm^2 -> m^2
    end
end