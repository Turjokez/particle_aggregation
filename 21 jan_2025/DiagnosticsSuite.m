classdef DiagnosticsTests
    %DIAGNOSTICSTESTS Numeric pass/fail checks for model bookkeeping.
    %
    % Usage:
    %   R = DiagnosticsTests.runAll();
    %
    % Requires rhs.decomposeTerms returning at least:
    %   dv_tot, dv_lin, dv_pp, dv_coag, dv_disagg
    % Updated CoagulationRHS adds:
    %   dv_other

    methods (Static)

        function R = runAll()
            fprintf('\n=== DiagnosticsTests: RUN ALL ===\n');

            R = struct();
            R.tstamp = string(datetime("now"));

            % ==========================================================
            % NEW (recommended): physically correct + fast tests
            % ==========================================================
            R.test_coag_conservation   = DiagnosticsTests.testCoagConservation_FAST_0D();
            R.test_disagg_conservation = DiagnosticsTests.testDisaggConservation_FAST_0D();
            R.test_linear_budget       = DiagnosticsTests.testLinearBudget_COLUMN();
            R.test_full_closure        = DiagnosticsTests.testFullClosure();
            R.test_no_hidden_coag      = DiagnosticsTests.testNoHiddenTerm_CoagOnly();

            % ==========================================================
            % OLD (kept): your original versions (optional)
            % Uncomment if you still want them:
            % ==========================================================
            % R.old_test_coag_conservation   = DiagnosticsTests.testCoagConservation();
            % R.old_test_disagg_conservation = DiagnosticsTests.testDisaggConservation();
            % R.old_test_linear_budget       = DiagnosticsTests.testLinearBudget();

            fprintf('\n=== DiagnosticsTests: SUMMARY ===\n');
            DiagnosticsTests.printResult("Coag conservation",   R.test_coag_conservation);
            DiagnosticsTests.printResult("Disagg conservation", R.test_disagg_conservation);
            DiagnosticsTests.printResult("Linear budget",       R.test_linear_budget);
            DiagnosticsTests.printResult("Full closure",        R.test_full_closure);
            DiagnosticsTests.printResult("No hidden (coag)",    R.test_no_hidden_coag);

            fprintf('=== DONE ===\n\n');
        end

        % ==========================================================
        % NEW TEST 1: Coag-only conservation in 0-D closed system
        % ==========================================================
        function out = testCoagConservation_FAST_0D()

            out = DiagnosticsTests.blankResult("testCoagConservation_FAST_0D");

            cfg = SimulationConfig();

            % Robustly force 0-D (different config versions)
            try, cfg.use_column = false; catch, end
            try, cfg.column = false; catch, end
            try, cfg.enable_column = false; catch, end

            % ---- isolate coag, kill everything else ----
            if isprop(cfg,'enable_coag'),    cfg.enable_coag   = true;  end
            if isprop(cfg,'enable_disagg'),  cfg.enable_disagg = false; end
            if isprop(cfg,'enable_sinking'), cfg.enable_sinking = false; end
            if isprop(cfg,'growth'),         cfg.growth        = 0; end
            if isprop(cfg,'enable_pp'),      cfg.enable_pp     = false; end

            % keep it fast
            if isprop(cfg,'t_max'),          cfg.t_max         = 2; end
            if isprop(cfg,'dt_out'),         cfg.dt_out        = 0.1; end

            sim = CoagulationSimulation(cfg);
            sim.run();

            [M, okM] = DiagnosticsTests.getInventory(sim);
            if ~okM
                out.ok = false;
                out.msg = "Inventory not available -> cannot test coag conservation";
                return;
            end

            M0   = M(1);
            Mend = M(end);

            rel = abs(Mend - M0) / max(abs(M0), 1e-30);

            out.metrics.M0 = M0;
            out.metrics.Mend = Mend;
            out.metrics.abs_dM = abs(Mend - M0);
            out.metrics.rel_dM = rel;

            tol_rel = 1e-3;
            out.ok = (rel < tol_rel);
            out.msg = sprintf("M0=%.3e Mend=%.3e | rel drift=%.3e (tol %.1e)", M0, Mend, rel, tol_rel);
        end

        % ==========================================================
        % NEW TEST 2: Disagg-only conservation in 0-D closed system
        % ==========================================================
        function out = testDisaggConservation_FAST_0D()

            out = DiagnosticsTests.blankResult("testDisaggConservation_FAST_0D");

            cfg = SimulationConfig();
            try, cfg.use_column = false; catch, end
            try, cfg.column = false; catch, end
            try, cfg.enable_column = false; catch, end

            % ---- isolate disagg ----
            if isprop(cfg,'enable_coag'),     cfg.enable_coag   = false; end
            if isprop(cfg,'enable_disagg'),   cfg.enable_disagg = true;  end
            if isprop(cfg,'disagg_apply_in'), cfg.disagg_apply_in = "rhs"; end
            if isprop(cfg,'disagg_rate'),     cfg.disagg_rate   = 1; end

            if isprop(cfg,'enable_sinking'),  cfg.enable_sinking = false; end
            if isprop(cfg,'growth'),          cfg.growth        = 0; end
            if isprop(cfg,'enable_pp'),       cfg.enable_pp     = false; end

            if isprop(cfg,'epsilon_const'),   cfg.epsilon_const = 1e-7; end
            if isprop(cfg,'eps_ref'),         cfg.eps_ref       = 1e-7; end

            if isprop(cfg,'t_max'),           cfg.t_max         = 2; end
            if isprop(cfg,'dt_out'),          cfg.dt_out        = 0.1; end

            sim = CoagulationSimulation(cfg);
            sim.run();

            [M, okM] = DiagnosticsTests.getInventory(sim);
            if ~okM
                out.ok = false;
                out.msg = "Inventory not available -> cannot test disagg conservation";
                return;
            end

            M0   = M(1);
            Mend = M(end);

            rel = abs(Mend - M0) / max(abs(M0), 1e-30);

            out.metrics.M0 = M0;
            out.metrics.Mend = Mend;
            out.metrics.abs_dM = abs(Mend - M0);
            out.metrics.rel_dM = rel;

            tol_rel = 1e-3;
            out.ok = (rel < tol_rel);
            out.msg = sprintf("M0=%.3e Mend=%.3e | rel drift=%.3e (tol %.1e)", M0, Mend, rel, tol_rel);
        end

        % ==========================================================
        % NEW TEST 3: COLUMN linear budget (FD vs integrated dv_lin)
        % ==========================================================
        function out = testLinearBudget_COLUMN()

            out = DiagnosticsTests.blankResult("testLinearBudget_COLUMN");

            cfg = SimulationConfig();
            try, cfg.use_column = true; catch, end
            try, cfg.column = true; catch, end
            try, cfg.enable_column = true; catch, end
            if isprop(cfg,'z_max'), cfg.z_max = 65; end
            if isprop(cfg,'dz'),    cfg.dz    = 5;  end

            if isprop(cfg,'enable_coag'),    cfg.enable_coag   = false; end
            if isprop(cfg,'enable_disagg'),  cfg.enable_disagg = false; end
            if isprop(cfg,'growth'),         cfg.growth        = 0; end
            if isprop(cfg,'enable_pp'),      cfg.enable_pp     = false; end

            if isprop(cfg,'t_max'),          cfg.t_max         = 5; end
            if isprop(cfg,'dt_out'),         cfg.dt_out        = 0.25; end

            sim = CoagulationSimulation(cfg);
            sim.run();

            t = sim.result.time(:);

            [M, okM] = DiagnosticsTests.getInventory(sim);
            if ~okM
                out.ok = false;
                out.msg = "Missing inventory -> cannot test linear budget";
                return;
            end

            dMdt_fd = gradient(M, t);

            have_terms = false;
            dMdt_lin = nan(size(t));

            try
                if ismethod(sim.rhs,'decomposeTerms')
                    have_terms = true;

                    cfg2 = sim.config;
                    Ns = cfg2.n_sections;
                    Nz = cfg2.getNumLayers();
                    dz_cm = cfg2.dz * 100;

                    for it = 1:numel(t)
                        v = sim.result.concentrations(it,:).';
                        terms = sim.rhs.decomposeTerms(t(it), v);
                        dMdt_lin(it) = DiagnosticsTests.integrateColumnRate(terms.dv_lin, Ns, Nz, dz_cm);
                    end
                end
            catch
                have_terms = false;
            end

            if have_terms && all(isfinite(dMdt_lin))
                resid = dMdt_fd - dMdt_lin;

                out.metrics.rms_resid = sqrt(mean(resid.^2));
                out.metrics.scale     = max(max(abs(dMdt_fd)), max(abs(dMdt_lin)));
                out.metrics.rel_rms   = out.metrics.rms_resid / max(out.metrics.scale, 1e-30);

                tol_rel = 2e-2;
                out.ok  = (out.metrics.rel_rms < tol_rel);
                out.msg = sprintf("FD vs ∫dv_lin dz: rel RMS resid = %.3e (tol %.1e)", out.metrics.rel_rms, tol_rel);
            else
                out.ok = false;
                out.msg = "decomposeTerms missing or dv_lin integration failed";
            end
        end

        % ==========================================================
        % TEST 4: Full closure INCLUDING dv_other
        % ==========================================================
        function out = testFullClosure()

            out = DiagnosticsTests.blankResult("testFullClosure");

            cfg = SimulationConfig();
            try, cfg.use_column = true; catch, end
            try, cfg.column = true; catch, end
            try, cfg.enable_column = true; catch, end
            if isprop(cfg,'z_max'), cfg.z_max = 65; end
            if isprop(cfg,'dz'),    cfg.dz    = 5;  end

            if isprop(cfg,'enable_coag'),      cfg.enable_coag   = true; end
            if isprop(cfg,'enable_disagg'),    cfg.enable_disagg = true; end
            if isprop(cfg,'disagg_apply_in'),  cfg.disagg_apply_in = "rhs"; end
            if isprop(cfg,'disagg_rate'),      cfg.disagg_rate   = 1; end
            if isprop(cfg,'epsilon_const'),    cfg.epsilon_const = 1e-7; end
            if isprop(cfg,'eps_ref'),          cfg.eps_ref       = 1e-7; end

            if isprop(cfg,'t_max'),            cfg.t_max         = 5; end
            if isprop(cfg,'dt_out'),           cfg.dt_out        = 0.25; end

            sim = CoagulationSimulation(cfg);
            sim.run();

            if ~ismethod(sim.rhs,'decomposeTerms')
                out.ok = false;
                out.msg = "rhs.decomposeTerms missing -> cannot test full closure";
                return;
            end

            t = sim.result.time(:);
            Y = sim.result.concentrations;

            ii = unique(round(linspace(1, numel(t), 6)));
            max_rel = 0;
            max_abs = 0;
            max_other_inf = 0;
            max_other_sum = 0;

            for k = 1:numel(ii)
                it = ii(k);
                v = Y(it,:).';

                terms = sim.rhs.decomposeTerms(t(it), v);

                if ~isfield(terms,'dv_other')
                    terms.dv_other = zeros(size(v)); % degrade gracefully
                end

                dv_sum = terms.dv_lin + terms.dv_coag + terms.dv_pp + terms.dv_disagg + terms.dv_other;
                dv_err = terms.dv_tot - dv_sum;

                abs_err = norm(dv_err, inf);
                scale   = max(1e-30, norm(terms.dv_tot, inf));
                rel_err = abs_err / scale;

                max_abs = max(max_abs, abs_err);
                max_rel = max(max_rel, rel_err);

                max_other_inf = max(max_other_inf, norm(terms.dv_other, inf));
                max_other_sum = max(max_other_sum, abs(sum(terms.dv_other)));
            end

            out.metrics.max_abs_inf = max_abs;
            out.metrics.max_rel_inf = max_rel;
            out.metrics.max_other_inf = max_other_inf;
            out.metrics.max_other_sum = max_other_sum;

            tol_rel = 1e-10;
            out.ok  = (max_rel < tol_rel);
            out.msg = sprintf("max rel inf error = %.3e (tol %.1e) | max||other||=%.3e | max|sum(other)|=%.3e", ...
                max_rel, tol_rel, max_other_inf, max_other_sum);
        end

        % ==========================================================
        % NEW TEST 5: No hidden term in coag-only 0-D
        % ==========================================================
        function out = testNoHiddenTerm_CoagOnly()

            out = DiagnosticsTests.blankResult("testNoHiddenTerm_CoagOnly");

            cfg = SimulationConfig();
            try, cfg.use_column = false; catch, end
            try, cfg.column = false; catch, end
            try, cfg.enable_column = false; catch, end

            if isprop(cfg,'enable_coag'),    cfg.enable_coag   = true;  end
            if isprop(cfg,'enable_disagg'),  cfg.enable_disagg = false; end
            if isprop(cfg,'enable_sinking'), cfg.enable_sinking = false; end
            if isprop(cfg,'growth'),         cfg.growth        = 0; end
            if isprop(cfg,'enable_pp'),      cfg.enable_pp     = false; end

            if isprop(cfg,'t_max'),          cfg.t_max         = 0.1; end
            if isprop(cfg,'dt_out'),         cfg.dt_out        = 0.1; end

            sim = CoagulationSimulation(cfg);
            sim.run();

            if ~ismethod(sim.rhs,'decomposeTerms')
                out.ok = false;
                out.msg = "rhs.decomposeTerms missing";
                return;
            end

            t0 = sim.result.time(1);
            v0 = sim.result.concentrations(1,:).';

            terms = sim.rhs.decomposeTerms(t0, v0);

            if ~isfield(terms,'dv_other')
                out.ok = false;
                out.msg = "dv_other missing (update CoagulationRHS.decomposeTerms first)";
                return;
            end

            other_sum = sum(terms.dv_other);
            other_inf = norm(terms.dv_other, inf);

            out.metrics.other_sum = other_sum;
            out.metrics.other_inf = other_inf;

            tol_sum = 1e-12;
            tol_inf = 1e-10;

            out.ok = (abs(other_sum) < tol_sum) && (other_inf < tol_inf);
            out.msg = sprintf("sum(other)=%.3e (tol %.1e) | inf(other)=%.3e (tol %.1e)", ...
                other_sum, tol_sum, other_inf, tol_inf);
        end

        % ==========================================================
        % OLD (KEPT) — your original tests (unchanged) BELOW
        % ==========================================================
        function out = testCoagConservation() %#ok<DEFNU>
            out = DiagnosticsTests.blankResult("testCoagConservation");
            cfg = SimulationConfig();
            cfg.use_column = true; cfg.z_max = 65; cfg.dz = 5;
            if isprop(cfg,'enable_coag'),  cfg.enable_coag  = true; end
            if isprop(cfg,'enable_disagg'),cfg.enable_disagg= false; end
            if isprop(cfg,'growth'),       cfg.growth       = 0; end
            if isprop(cfg,'enable_pp'),    cfg.enable_pp    = false; end
            sim = CoagulationSimulation(cfg);
            sim.run();
            [M, okM] = DiagnosticsTests.getInventory(sim);
            if ~okM
                out.ok = false;
                out.msg = "Inventory not available -> cannot test conservation";
                return;
            end
            dM = max(M) - min(M);
            rel = dM / max(max(M), 1e-30);
            out.metrics.abs_dM = dM;
            out.metrics.rel_dM = rel;
            tol_rel = 1e-3;
            out.ok = (rel < tol_rel);
            out.msg = sprintf("rel drift = %.3e (tol %.1e)", rel, tol_rel);
        end

        function out = testDisaggConservation() %#ok<DEFNU>
            out = DiagnosticsTests.blankResult("testDisaggConservation");
            cfg = SimulationConfig();
            cfg.use_column = true; cfg.z_max = 65; cfg.dz = 5;
            if isprop(cfg,'enable_coag'),   cfg.enable_coag   = false; end
            if isprop(cfg,'enable_disagg'), cfg.enable_disagg = true; end
            if isprop(cfg,'disagg_apply_in'), cfg.disagg_apply_in = "rhs"; end
            if isprop(cfg,'disagg_rate'),   cfg.disagg_rate   = 1; end
            if isprop(cfg,'growth'),        cfg.growth        = 0; end
            if isprop(cfg,'enable_pp'),     cfg.enable_pp     = false; end
            if isprop(cfg,'epsilon_const'), cfg.epsilon_const = 1e-7; end
            if isprop(cfg,'eps_ref'),       cfg.eps_ref       = 1e-7; end
            sim = CoagulationSimulation(cfg);
            sim.run();
            [M, okM] = DiagnosticsTests.getInventory(sim);
            if ~okM
                out.ok = false;
                out.msg = "Inventory not available -> cannot test conservation";
                return;
            end
            dM = max(M) - min(M);
            rel = dM / max(max(M), 1e-30);
            out.metrics.abs_dM = dM;
            out.metrics.rel_dM = rel;
            tol_rel = 1e-3;
            out.ok = (rel < tol_rel);
            out.msg = sprintf("rel drift = %.3e (tol %.1e)", rel, tol_rel);
        end

        function out = testLinearBudget() %#ok<DEFNU>
            out = DiagnosticsTests.blankResult("testLinearBudget");
            cfg = SimulationConfig();
            cfg.use_column = true; cfg.z_max = 65; cfg.dz = 5;
            if isprop(cfg,'enable_coag'),   cfg.enable_coag   = false; end
            if isprop(cfg,'enable_disagg'), cfg.enable_disagg = false; end
            if isprop(cfg,'growth'),        cfg.growth        = 0; end
            if isprop(cfg,'enable_pp'),     cfg.enable_pp     = false; end
            sim = CoagulationSimulation(cfg);
            sim.run();
            t = sim.result.time(:);
            [M, okM] = DiagnosticsTests.getInventory(sim);
            [F, okF] = DiagnosticsTests.getBottomExport(sim);
            if ~(okM && okF)
                out.ok = false;
                out.msg = "Missing inventory/export -> cannot test linear budget";
                return;
            end
            dMdt = gradient(M, t);
            resid = dMdt + F;
            out.metrics.max_abs_resid = max(abs(resid));
            out.metrics.rms_resid     = sqrt(mean(resid.^2));
            out.metrics.scale         = max(max(abs(dMdt)), max(abs(F)));
            out.metrics.rel_rms       = out.metrics.rms_resid / max(out.metrics.scale, 1e-30);
            tol_rel = 2e-2;
            out.ok  = (out.metrics.rel_rms < tol_rel);
            out.msg = sprintf("rel RMS resid = %.3e (tol %.1e)", out.metrics.rel_rms, tol_rel);
        end

        % ==========================================================
        % Helpers
        % ==========================================================
        function [M, ok] = getInventory(sim)
            ok = false; M = [];
            try
                out = sim.result;
                grd = sim.grid;
                cfg = sim.config;

                if ~isfield(out,'output_data') || isempty(out.output_data)
                    out.output_data = OutputGenerator.spectraAndFluxes(out.time, out.concentrations, grd, cfg);
                end
                od = out.output_data;

                if isfield(od,'column_total_inventory_cm3m2') && ~isempty(od.column_total_inventory_cm3m2)
                    M = od.column_total_inventory_cm3m2(:);
                    ok = true;
                    return;
                end
                if isfield(od,'total_mass') && ~isempty(od.total_mass)
                    M = od.total_mass(:);
                    ok = true;
                    return;
                end
            catch
            end
        end

        function [F, ok] = getBottomExport(sim)
            ok = false; F = [];
            try
                out = sim.result;
                grd = sim.grid;
                cfg = sim.config;

                if ~isfield(out,'output_data') || isempty(out.output_data)
                    out.output_data = OutputGenerator.spectraAndFluxes(out.time, out.concentrations, grd, cfg);
                end
                od = out.output_data;

                if isfield(od,'bottom_total_flux_cm3m2d') && ~isempty(od.bottom_total_flux_cm3m2d)
                    F = od.bottom_total_flux_cm3m2d(:);
                    ok = true;
                    return;
                end
                if isfield(od,'total_flux') && ~isempty(od.total_flux)
                    F = od.total_flux(:);
                    ok = true;
                    return;
                end
            catch
            end
        end

        function rate_cm3m2d = integrateColumnRate(dv, Ns, Nz, dz_cm)
            N2 = reshape(dv(:), [Ns, Nz]);
            colsum = sum(N2, 1);
            rate_cm3m2d = sum(colsum) * dz_cm * 1e4;
        end

        function out = blankResult(name)
            out = struct();
            out.name = string(name);
            out.ok = false;
            out.msg = "";
            out.metrics = struct();
        end

        function printResult(label, r)
            if r.ok
                fprintf('[PASS] %-20s | %s\n', label, r.msg);
            else
                fprintf('[FAIL] %-20s | %s\n', label, r.msg);
            end
        end

    end
end