classdef DiagnosticsTests
    %DIAGNOSTICSTESTS Numeric pass/fail checks for model bookkeeping.
    %
    % Usage:
    %   clear classes;
    %   R = DiagnosticsTests.runAll();
    %
    % Assumptions:
    % - sim.run() populates sim.result with .time, .concentrations
    % - RHS implements decomposeTerms returning struct:
    %     dv_tot, dv_lin, dv_coag, dv_pp, dv_disagg (+ optional dv_other)
    %
    % Key point:
    % - Conservation depends on what the STATE means.
    %   We compute and report THREE weights:
    %       (0) ones  : tests conservation of sum(v)
    %       (1) Vbin  : tests conservation of sum(v .* Vbin)
    %       (2) av_vol: diagnostic only
    %
    % PASS criterion used here
    % - ones (because coag conserves ones, not Vbin)

    methods (Static)

        % ==========================================================
        % ENTRY POINT
        % ==========================================================
        function R = runAll()
            fprintf('\n=== DiagnosticsTests: RUN ALL ===\n');

            R = struct();
            R.tstamp = string(datetime("now"));

            R.test_coag_conservation   = DiagnosticsTests.testCoagConservation_FAST_0D();
            R.test_disagg_conservation = DiagnosticsTests.testDisaggConservation_FAST_0D();
            R.test_linear_budget       = DiagnosticsTests.testLinearBudget_COLUMN();
            R.test_full_closure        = DiagnosticsTests.testFullClosure();
            R.test_no_hidden_coag      = DiagnosticsTests.testNoHiddenTerm_CoagOnly();

            fprintf('\n=== DiagnosticsTests: SUMMARY ===\n');
            DiagnosticsTests.printResult("Coag conservation",   R.test_coag_conservation);
            DiagnosticsTests.printResult("Disagg conservation", R.test_disagg_conservation);
            DiagnosticsTests.printResult("Linear budget",       R.test_linear_budget);
            DiagnosticsTests.printResult("Full closure",        R.test_full_closure);
            DiagnosticsTests.printResult("No hidden (coag)",    R.test_no_hidden_coag);

            fprintf('=== DONE ===\n\n');
        end

        % ==========================================================
        % TEST 1: Coag-only conservation (0-D)
        % ==========================================================
        function out = testCoagConservation_FAST_0D()

            out = DiagnosticsTests.blankResult("testCoagConservation_FAST_0D");

            cfg = SimulationConfig();
            cfg.use_column = false;

            cfg.enable_coag    = true;
            cfg.enable_disagg  = false;
            cfg.enable_sinking = false;
            cfg.enable_pp      = false;
            cfg.growth         = 0;

            cfg.t_init  = 0;
            cfg.t_final = 0.2;
            cfg.delta_t = 0.2;

            sim = CoagulationSimulation(cfg);
            sim.run();

            if ~ismethod(sim.rhs,'decomposeTerms')
                out.ok  = false;
                out.msg = "rhs.decomposeTerms missing";
                return;
            end

            t0 = sim.result.time(1);
            v0 = sim.result.concentrations(1,:).';
            terms = sim.rhs.decomposeTerms(t0, v0);

            if ~isfield(terms,'dv_coag') || isempty(terms.dv_coag)
                out.ok  = false;
                out.msg = "dv_coag missing/empty";
                return;
            end

            W = DiagnosticsTests.getWeightsForConservation(sim);

            Ns = cfg.n_sections;
            Nz = 1;
            dz_cm = cfg.dz * 100;

            out.metrics.ones_relRate_per_day  = DiagnosticsTests.relRate(v0, terms.dv_coag, W, "ones",  Ns, Nz, dz_cm);
            out.metrics.vbin_relRate_per_day  = DiagnosticsTests.relRate(v0, terms.dv_coag, W, "vbin",  Ns, Nz, dz_cm);
            out.metrics.avvol_relRate_per_day = DiagnosticsTests.relRate(v0, terms.dv_coag, W, "avvol", Ns, Nz, dz_cm);

            tol = 3e-4;
            out.ok = isfinite(out.metrics.ones_relRate_per_day) && abs(out.metrics.ones_relRate_per_day) < tol;

            out.msg = sprintf("PASS criterion: ones | ones relRate=%.3e /d (tol %.1e)", ...
                              out.metrics.ones_relRate_per_day, tol);
        end

        % ==========================================================
        % TEST 2: Disagg-only conservation (0-D)
        % ==========================================================
        function out = testDisaggConservation_FAST_0D()

            out = DiagnosticsTests.blankResult("testDisaggConservation_FAST_0D");

            cfg = SimulationConfig();
            cfg.use_column = false;

            cfg.enable_coag     = false;
            cfg.enable_disagg   = true;
            cfg.disagg_apply_in = "rhs";
            cfg.disagg_rate     = 1;
            cfg.enable_pp       = false;
            cfg.growth          = 0;

            cfg.epsilon_const = 1e-7;
            cfg.eps_ref       = 1e-7;

            cfg.t_init  = 0;
            cfg.t_final = 0.2;
            cfg.delta_t = 0.2;

            sim = CoagulationSimulation(cfg);
            sim.run();

            if ~ismethod(sim.rhs,'decomposeTerms')
                out.ok  = false;
                out.msg = "rhs.decomposeTerms missing";
                return;
            end

            t0 = sim.result.time(1);
            v0 = sim.result.concentrations(1,:).';
            terms = sim.rhs.decomposeTerms(t0, v0);

            if ~isfield(terms,'dv_disagg') || isempty(terms.dv_disagg)
                out.ok  = false;
                out.msg = "dv_disagg missing/empty";
                return;
            end

            W = DiagnosticsTests.getWeightsForConservation(sim);

            Ns = cfg.n_sections;
            Nz = 1;
            dz_cm = cfg.dz * 100;

            out.metrics.ones_relRate_per_day = DiagnosticsTests.relRate(v0, terms.dv_disagg, W, "ones", Ns, Nz, dz_cm);

            tol = 1e-10;
            out.ok = isfinite(out.metrics.ones_relRate_per_day) && abs(out.metrics.ones_relRate_per_day) < tol;

            out.msg = sprintf("PASS criterion: ones | ones relRate=%.3e /d (tol %.1e)", ...
                              out.metrics.ones_relRate_per_day, tol);
        end

        % ==========================================================
        % TEST 3: Column linear budget (FD vs ∫dv_lin)
        % ==========================================================
        function out = testLinearBudget_COLUMN()

            out = DiagnosticsTests.blankResult("testLinearBudget_COLUMN");

            cfg = SimulationConfig();
            cfg.use_column = true;
            cfg.z_max = 65;
            cfg.dz    = 5;

            cfg.enable_coag   = false;
            cfg.enable_disagg = false;
            cfg.enable_pp     = false;
            cfg.growth        = 0;

            cfg.t_init  = 0;
            cfg.t_final = 5;
            cfg.delta_t = 0.25;

            sim = CoagulationSimulation(cfg);
            sim.run();

            if ~ismethod(sim.rhs,'decomposeTerms')
                out.ok  = false;
                out.msg = "rhs.decomposeTerms missing";
                return;
            end

            t = sim.result.time(:);
            Y = sim.result.concentrations;

            W = DiagnosticsTests.getWeightsForConservation(sim);
            w = W.wOnes;

            Ns = cfg.n_sections;
            Nz = cfg.getNumLayers();
            dz_cm = cfg.dz * 100;

            nt = numel(t);
            M = zeros(nt,1);
            dMdt_lin = zeros(nt,1);

            for it = 1:nt
                v = Y(it,:).';
                M(it) = DiagnosticsTests.integrateState(v, w, Ns, Nz, dz_cm);

                T = sim.rhs.decomposeTerms(t(it), v);
                if ~isfield(T,'dv_lin') || isempty(T.dv_lin)
                    out.ok  = false;
                    out.msg = "dv_lin missing/empty";
                    return;
                end
                dMdt_lin(it) = DiagnosticsTests.integrateState(T.dv_lin, w, Ns, Nz, dz_cm);
            end

            % central difference FD
            dMdt_fd = nan(nt,1);
            for i = 2:nt-1
                dMdt_fd(i) = (M(i+1) - M(i-1)) / (t(i+1) - t(i-1));
            end

            ok = isfinite(dMdt_fd) & isfinite(dMdt_lin);
            ok([1 end]) = false;

            resid = dMdt_fd(ok) - dMdt_lin(ok);

            rms_resid = sqrt(mean(resid.^2));
            scale = max([max(abs(dMdt_fd(ok))) max(abs(dMdt_lin(ok)))]);
            rel_rms = rms_resid / max(scale,1e-30);

            tol = 2e-2;
            out.ok = isfinite(rel_rms) && rel_rms < tol;
            out.msg = sprintf("FD vs ∫dv_lin (ones): rel RMS resid = %.3e (tol %.1e)", rel_rms, tol);
        end

        % ==========================================================
        % TEST 4: Full closure
        % ==========================================================
        function out = testFullClosure()

            out = DiagnosticsTests.blankResult("testFullClosure");

            cfg = SimulationConfig();
            cfg.use_column = true;
            cfg.z_max = 65;
            cfg.dz    = 5;

            cfg.enable_coag     = true;
            cfg.enable_disagg   = true;
            cfg.disagg_apply_in = "rhs";
            cfg.disagg_rate     = 1;
            cfg.epsilon_const   = 1e-7;
            cfg.eps_ref         = 1e-7;

            cfg.t_init  = 0;
            cfg.t_final = 5;
            cfg.delta_t = 0.25;

            sim = CoagulationSimulation(cfg);
            sim.run();

            if ~ismethod(sim.rhs,'decomposeTerms')
                out.ok  = false;
                out.msg = "rhs.decomposeTerms missing";
                return;
            end

            t = sim.result.time(:);
            Y = sim.result.concentrations;

            ii = unique(round(linspace(1, numel(t), 6)));
            max_rel = 0;

            for k = ii
                v = Y(k,:).';
                T = sim.rhs.decomposeTerms(t(k), v);
                if ~isfield(T,'dv_other') || isempty(T.dv_other)
                    T.dv_other = zeros(size(v));
                end

                % Guard missing fields
                if ~isfield(T,'dv_lin'),    T.dv_lin    = zeros(size(v)); end
                if ~isfield(T,'dv_coag'),   T.dv_coag   = zeros(size(v)); end
                if ~isfield(T,'dv_pp'),     T.dv_pp     = zeros(size(v)); end
                if ~isfield(T,'dv_disagg'), T.dv_disagg = zeros(size(v)); end

                dv_sum = T.dv_lin + T.dv_coag + T.dv_pp + T.dv_disagg + T.dv_other;
                dv_err = T.dv_tot - dv_sum;

                rel = norm(dv_err,inf) / max(norm(T.dv_tot,inf),1e-30);
                max_rel = max(max_rel, rel);
            end

            tol = 1e-10;
            out.ok = isfinite(max_rel) && max_rel < tol;
            out.msg = sprintf("max rel inf error = %.3e (tol %.1e)", max_rel, tol);
        end

        % ==========================================================
        % TEST 5: No hidden term in coag-only
        % ==========================================================
        function out = testNoHiddenTerm_CoagOnly()

            out = DiagnosticsTests.blankResult("testNoHiddenTerm_CoagOnly");

            cfg = SimulationConfig();
            cfg.use_column = false;

            cfg.enable_coag   = true;
            cfg.enable_disagg = false;
            cfg.enable_pp     = false;
            cfg.growth        = 0;

            cfg.t_init  = 0;
            cfg.t_final = 0.1;
            cfg.delta_t = 0.1;

            sim = CoagulationSimulation(cfg);
            sim.run();

            if ~ismethod(sim.rhs,'decomposeTerms')
                out.ok  = false;
                out.msg = "rhs.decomposeTerms missing";
                return;
            end

            v0 = sim.result.concentrations(1,:).';
            T  = sim.rhs.decomposeTerms(sim.result.time(1), v0);

            if ~isfield(T,'dv_other') || isempty(T.dv_other)
                T.dv_other = zeros(size(v0));
            end

            tol_sum = 1e-12;
            tol_inf = 1e-10;

            out.ok = abs(sum(T.dv_other)) < tol_sum && norm(T.dv_other,inf) < tol_inf;
            out.msg = sprintf("sum(other)=%.3e | inf(other)=%.3e", sum(T.dv_other), norm(T.dv_other,inf));
        end

        % ==========================================================
        % HELPERS
        % ==========================================================
        function rr = relRate(v0, dv, W, which, Ns, Nz, dz_cm)
            rr = nan;
            try
                if strcmpi(which,"ones")
                    w = W.wOnes;
                elseif strcmpi(which,"vbin")
                    w = W.wVbin;
                elseif strcmpi(which,"avvol")
                    w = W.wAvVol;
                else
                    return;
                end
                T0 = DiagnosticsTests.integrateState(v0, w, Ns, Nz, dz_cm);
                dT = DiagnosticsTests.integrateState(dv, w, Ns, Nz, dz_cm);
                rr = dT / max(abs(T0),1e-30);
            catch
            end
        end

        function I = integrateState(vflat, w, Ns, Nz, dz_cm)
            if Nz==1
                I = sum(vflat(:).*w(:));
            else
                N2 = reshape(vflat,[Ns Nz]);
                W2 = repmat(w(:),1,Nz);
                I = sum(sum(N2.*W2))*dz_cm;
            end
        end

        function W = getWeightsForConservation(sim)
            % Robust weight builder:
            % - always provide ones
            % - Vbin: try grid method or Disaggregation helper, else sphere fallback
            % - av_vol: try sim.grid.av_vol, else fallback to Vbin
            Ns = sim.config.n_sections;

            W = struct();
            W.wOnes = ones(Ns,1);

            % --- Vbin ---
            wV = [];
            try
                gridS = sim.rhs.getSectionGrid();
                if ismethod(gridS,'getBinVolumes')
                    wV = gridS.getBinVolumes();
                    wV = wV(:);
                end
            catch
            end
            if isempty(wV)
                try
                    wV = Disaggregation.getBinVolumes_cm3(sim.rhs.getSectionGrid());
                    wV = wV(:);
                catch
                end
            end
            if isempty(wV)
                % fallback: just use ones (or sphere if you want)
                wV = ones(Ns,1);
            end
            W.wVbin = wV;

            % --- av_vol ---
            wA = [];
            try
                if isprop(sim.grid,'av_vol') && ~isempty(sim.grid.av_vol)
                    wA = sim.grid.av_vol(:);
                end
            catch
            end
            if isempty(wA)
                wA = W.wVbin;
            end
            W.wAvVol = wA;

            W.ok_any = true;
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