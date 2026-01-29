classdef CoagulationRHS < handle
    % COAGULATIONRHS ODE right-hand side for coagulation equations
    %
    % coagulation + growth + sinking + disaggregation.
    %
    % NEW-2025-12-21 (CRITICAL UNITS FIX):
    % Solver time is in DAYS, so RHS must return dv/dt in [per day].
    %
    % NEW-2026-01-11 (DISAGG FIX):
    % Disaggregation.applyWithScaling is an instantaneous mapping Y -> Ynew.
    % To use inside an ODE RHS, convert into a tendency:
    %   dY/dt|disagg = disagg_rate * (Y_after - Y_before)   [per day]
    %
    % NEW-2026-01-11 (GRID FIX):
    % Adds helper getSectionGrid() and use it everywhere we call
    % Disaggregation.applyWithScaling.
    %
    % NEW-2026-01-14 (BETA ORIENTATION SWITCH):
    % Optional: transpose ONLY b1 and b3 in RHS evaluation.
    % Controlled by:
    %   cfg.beta_row_form = "row_nrm"   (preferred)
    % OR
    %   cfg.beta_fix_b1_b3 = true
    %
    % NEW-2026-01-14 (STEP-1 BUGFIX):
    % Legacy disagg matrices (disaggMinus/disaggPlus) must NOT act when
    % cfg.enable_disagg == false.
    %
    % NEW-2026-01-15 (PP SOURCE FIXED FOR AREAL UNITS):
    % Adds explicit PP source injection dv_pp into (pp_bin, pp_layer).
    %
    % IMPORTANT BEHAVIOR (NO pp_rate_units PROPERTY REQUIRED):
    % - If cfg.export_weight == "vbin": cfg.pp_rate is interpreted as
    %       pp_rate_areal = [cm^3 / cm^2 / day]
    %   and we convert to a concentration tendency (state/day) as:
    %       dv_pp = pp_rate_areal / (wbin(pp_bin) * dz_cm)
    %   where wbin is bin volume [cm^3/particle], dz_cm is layer thickness [cm].
    %
    % - Otherwise: cfg.pp_rate is interpreted as already being in state units/day.
    %
    % NEW-2026-01-XX (COAG BV CLOSURE):
    % If cfg.enforce_coag_bv_closure == true, enforce:
    %   sum_s ( Vbin(s) * dv_coag(s) ) = 0
    % per layer (and in 0-D), by applying a minimal correction to the last bin.
    % This fixes leakage caused by overflow pairs (Vp > Vmax) that get clipped.
    %
    % NOTE:
    % This file supports BOTH:
    %   terms = rhs.decomposeTerms(t,v)          % struct
    %   [dv_tot,dv_lin,dv_coag,dv_pp,dv_disagg,dv_other] = rhs.decomposeTerms(t,v)

    properties
        betas;          % BetaMatrices object
        linear;         % Linear operator (0-D) OR column operator (Ns*Nz x Ns*Nz)
        disaggMinus;    % Disaggregation loss matrix (Ns x Ns)
        disaggPlus;     % Disaggregation gain matrix (Ns x Ns)
        config;         % SimulationConfig object

        eps_time = [];     % time vector [d]
        eps_vals = [];     % epsilon(t) [vector OR Nt x Nz]  (stored as Nt x Nz when possible)
        eps_const = [];    % fallback constant epsilon

        grid = [];         % DerivedGrid (optional; REQUIRED for vbin PP conversion and new disagg)
        z_cache = [];      % cached z-centers (m)
    end

    methods
        function obj = CoagulationRHS(betas, linear, disaggMinus, disaggPlus, config, varargin)
            obj.betas       = betas;
            obj.linear      = linear;
            obj.disaggMinus = disaggMinus;
            obj.disaggPlus  = disaggPlus;
            obj.config      = config;

            if isprop(config, 'epsilon_const')
                obj.eps_const = config.epsilon_const;
            else
                obj.eps_const = [];
            end

            % allow passing grid as 6th input (no breaking older calls)
            if nargin >= 6 && ~isempty(varargin)
                try
                    obj.grid = varargin{1};
                catch
                    obj.grid = [];
                end
            end

            % cache depth grid if possible
            try
                if ismethod(config,'getZ')
                    obj.z_cache = config.getZ();
                end
            catch
                obj.z_cache = [];
            end
        end

        % ==========================================================
        % NEW-2026-01-14: beta orientation helper (safe default)
        % ==========================================================
        function [b1,b2,b3,b4,b5] = getBetasForRHS(obj)
            b1 = obj.betas.b1;
            b2 = obj.betas.b2;
            b3 = obj.betas.b3;
            b4 = obj.betas.b4;
            b5 = obj.betas.b5;

            mode = "legacy";

            try
                if isprop(obj.config,'beta_row_form') && ~isempty(obj.config.beta_row_form)
                    mode = lower(string(obj.config.beta_row_form));
                end
            catch
            end

            try
                if mode == "legacy"
                    if isprop(obj.config,'beta_fix_b1_b3') && ~isempty(obj.config.beta_fix_b1_b3) && logical(obj.config.beta_fix_b1_b3)
                        mode = "row_nrm";
                    end
                end
            catch
            end

            if mode == "row_nrm" || mode == "row_nr_m" || mode == "row"
                b1 = b1.';   % transpose b1
                b3 = b3.';   % transpose b3
            end
        end

        % ==========================================================
        % NEW-2026-01-XX: coag biovolume closure helper
        % ==========================================================
        function dvk = applyCoagBvClosure(obj, dvk)
            do_closure = false;
            try
                if isprop(obj.config,'enforce_coag_bv_closure') && ~isempty(obj.config.enforce_coag_bv_closure)
                    do_closure = logical(obj.config.enforce_coag_bv_closure);
                end
            catch
                do_closure = false;
            end
            if ~do_closure
                return;
            end

            V = [];
            try
                if ~isempty(obj.grid) && isprop(obj.grid,'av_vol') && ~isempty(obj.grid.av_vol)
                    V = obj.grid.av_vol(:); % cm^3 per particle
                end
            catch
                V = [];
            end
            if isempty(V)
                % do not crash; just skip
                return;
            end

            dvk = dvk(:);
            if numel(dvk) ~= numel(V)
                return;
            end

            leak = sum(dvk .* V); % (state/day)*(cm^3/particle)
            if ~isfinite(leak) || leak == 0
                return;
            end

            if V(end) > 0
                dvk(end) = dvk(end) - leak / V(end);
            end
        end

        function setEpsilonSeries(obj, t, eps_vals)
            obj.eps_time = t(:);

            if isvector(eps_vals)
                obj.eps_vals = eps_vals(:);
                Nt = numel(obj.eps_time);
                Ne = numel(obj.eps_vals);
                n  = min(Nt, Ne);
                obj.eps_time = obj.eps_time(1:n);
                obj.eps_vals = obj.eps_vals(1:n);
                return;
            end

            Nt = numel(obj.eps_time);

            % Case A: already Nt x Nz
            if size(eps_vals,1) == Nt
                obj.eps_vals = eps_vals;
                return;
            end

            % Case B: Nz x Nt -> transpose to Nt x Nz
            if size(eps_vals,2) == Nt && size(eps_vals,1) ~= Nt
                obj.eps_vals = eps_vals.';
                return;
            end

            obj.eps_vals = eps_vals; % fallback
        end

        function setGrid(obj, gridObj)
            obj.grid = gridObj;
        end

        % --------------------------------------------------------------
        % Always return a valid "section grid" for Disaggregation
        % --------------------------------------------------------------
        function gridS = getSectionGrid(obj)
            gridS = obj.grid;

            if isempty(gridS)
                error('CoagulationRHS:getSectionGrid:EmptyGrid', ...
                    'obj.grid is empty. Pass DerivedGrid into CoagulationRHS constructor or call rhs.setGrid(grid).');
            end

            try
                if ismethod(gridS,'getConservedRadii') || ismethod(gridS,'getFractalRadii')
                    return;
                end
            catch
            end

            rv = [];
            try
                if isprop(gridS,'av_vol') && ~isempty(gridS.av_vol)
                    av = gridS.av_vol(:);
                    rv = (0.75/pi * av).^(1/3);
                end
            catch
            end

            if isempty(rv)
                error('CoagulationRHS:getSectionGrid:InvalidGrid', ...
                    'obj.grid does not implement getConservedRadii/getFractalRadii and no fallback radii found.');
            end

            gridS = SectionGridAdapter(rv, obj.config);
        end

        function dvdt = rhs(obj, t, v)
            dvdt = obj.evaluate(t, v);

            if any(~isfinite(dvdt))
                error('CoagulationRHS:NaNInf', ...
                      'RHS produced NaN/Inf at t=%.6g. min=%g max=%g', ...
                      t, min(dvdt(~isnan(dvdt))), max(dvdt(~isnan(dvdt))));
            end
        end

        function eps_now = getEpsScalar(obj, t)
            eps_now = [];

            if ~isempty(obj.eps_time) && ~isempty(obj.eps_vals)

                if isvector(obj.eps_vals)
                    Nt = numel(obj.eps_time);
                    Ne = numel(obj.eps_vals);
                    n  = min(Nt, Ne);
                    tt = obj.eps_time(1:n);
                    ee = obj.eps_vals(1:n);
                    try
                        eps_now = interp1(tt(:), ee(:), t, 'linear', 'extrap');
                    catch
                        eps_now = [];
                    end

                else
                    try
                        tt = obj.eps_time(:);
                        Em = obj.eps_vals;

                        if size(Em,1) == numel(tt)          % Nt x Nz
                            eps_now = interp1(tt, Em(:,1), t, 'linear', 'extrap');
                        elseif size(Em,2) == numel(tt)      % Nz x Nt
                            eps_now = interp1(tt, Em(1,:).', t, 'linear', 'extrap');
                        end
                    catch
                        eps_now = [];
                    end
                end
            end

            if isempty(eps_now) && ~isempty(obj.eps_const)
                eps_now = obj.eps_const;
            end

            try
                if ~isempty(eps_now) && isnumeric(eps_now) && isfinite(eps_now)
                    eps_now = max(eps_now, 0);
                end
            catch
            end
        end

        % ==========================================================
        % Epsilon at a specific layer (robust Nz x Nt OR Nt x Nz)
        % ==========================================================
        function eps_here = getEpsAtLayer(obj, t, k, zc_k)
            eps_here = [];

            % 1) eps_fun(t,z)
            try
                if isprop(obj.config,'eps_fun') && ~isempty(obj.config.eps_fun) && isa(obj.config.eps_fun,'function_handle')
                    eps_here = obj.config.eps_fun(t, zc_k);
                end
            catch
                eps_here = [];
            end

            if ~isempty(eps_here) && isfinite(eps_here)
                eps_here = max(eps_here, 0);
                return;
            end

            % 2) epsilon_time + epsilon_series
            try
                if isprop(obj.config,'epsilon_time') && isprop(obj.config,'epsilon_series') ...
                        && ~isempty(obj.config.epsilon_time) && ~isempty(obj.config.epsilon_series)

                    tt = obj.config.epsilon_time(:);
                    Em = obj.config.epsilon_series;

                    if ~isvector(Em) && size(Em,1) == numel(tt)       % Nt x Nz
                        if k <= size(Em,2)
                            eps_here = interp1(tt, Em(:,k), t, 'linear', 'extrap');
                        end
                    elseif ~isvector(Em) && size(Em,2) == numel(tt)   % Nz x Nt
                        if k <= size(Em,1)
                            eps_here = interp1(tt, Em(k,:).', t, 'linear', 'extrap');
                        end
                    end
                end
            catch
                eps_here = [];
            end

            if ~isempty(eps_here) && isfinite(eps_here)
                eps_here = max(eps_here, 0);
                return;
            end

            % 3) fallback scalar
            eps_here = obj.getEpsScalar(t);
            if ~isempty(eps_here) && isfinite(eps_here)
                eps_here = max(eps_here, 0);
            end
        end

        % ==========================================================
        % PP source helper (FIXED)
        % ==========================================================
        function dv_pp = getPPSourceVector(obj, v)
            v = v(:);
            dv_pp = zeros(size(v));

            enable_pp = false;
            try
                if isprop(obj.config,'enable_pp') && ~isempty(obj.config.enable_pp)
                    enable_pp = logical(obj.config.enable_pp);
                end
            catch
                enable_pp = false;
            end
            if ~enable_pp
                return;
            end

            % rate (meaning depends on export_weight)
            pp_rate = 0.0;
            try
                if isprop(obj.config,'pp_rate') && ~isempty(obj.config.pp_rate)
                    pp_rate = obj.config.pp_rate;
                end
            catch
                pp_rate = 0.0;
            end
            if ~isfinite(pp_rate) || pp_rate <= 0
                return;
            end

            % indices
            Ns = obj.config.n_sections;

            pp_bin = 1;
            try
                if isprop(obj.config,'pp_bin') && ~isempty(obj.config.pp_bin)
                    pp_bin = round(obj.config.pp_bin);
                end
            catch
                pp_bin = 1;
            end
            pp_bin = max(1, min(Ns, pp_bin));

            % decide interpretation of pp_rate
            is_vbin = false;
            try
                if isprop(obj.config,'export_weight') && ~isempty(obj.config.export_weight)
                    is_vbin = strcmpi(string(obj.config.export_weight),"vbin");
                end
            catch
                is_vbin = false;
            end

            if isprop(obj.config,'use_column') && obj.config.use_column
                Nz = obj.config.getNumLayers();

                if numel(v) ~= Ns * Nz
                    error('PP source: state length mismatch. Expected %d, got %d.', Ns*Nz, numel(v));
                end

                pp_layer = 1;
                try
                    if isprop(obj.config,'pp_layer') && ~isempty(obj.config.pp_layer)
                        pp_layer = round(obj.config.pp_layer);
                    end
                catch
                    pp_layer = 1;
                end
                pp_layer = max(1, min(Nz, pp_layer));

                idx = (pp_layer-1)*Ns + pp_bin;

                if is_vbin
                    % pp_rate is AREAL biovolume production: [cm^3/cm^2/day]
                    if isempty(obj.grid)
                        error('PP source (vbin): obj.grid is empty. Pass grid into CoagulationRHS or call rhs.setGrid(grid).');
                    end

                    dz_cm = obj.config.dz * 100;

                    % wbin (cm^3/particle)
                    try
                        wbin = obj.grid.getBinVolumes();
                    catch
                        wbin = Disaggregation.getBinVolumes_cm3(obj.grid);
                    end
                    wbin = wbin(:);

                    dv_pp(idx) = dv_pp(idx) + pp_rate / (wbin(pp_bin) * dz_cm);  % state/day
                else
                    % pp_rate is already in state units/day
                    dv_pp(idx) = dv_pp(idx) + pp_rate;
                end

            else
                % 0-D
                if is_vbin
                    error('PP source (vbin): 0-D mode not supported for areal pp_rate.');
                end
                if pp_bin >= 1 && pp_bin <= numel(v)
                    dv_pp(pp_bin) = dv_pp(pp_bin) + pp_rate;
                end
            end
        end

        function dvdt = evaluate(obj, t, v)
            Ns = obj.config.n_sections;

            % coag switch
            enable_coag = true;
            if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag)
                enable_coag = logical(obj.config.enable_coag);
            end

            % disagg switch (legacy matrices must obey enable_disagg)
            enable_disagg = false;
            try
                if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                    enable_disagg = logical(obj.config.enable_disagg);
                end
            catch
                enable_disagg = false;
            end

            % epsilon(t) (optional)
            eps_now = obj.getEpsScalar(t); %#ok<NASGU>

            % negative handling
            clip_negative = true;
            if isprop(obj.config, 'clip_negative')
                clip_negative = logical(obj.config.clip_negative);
            end

            % new disagg switch
            use_new_disagg = false;
            if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                if logical(obj.config.enable_disagg)
                    if isprop(obj.config,'disagg_apply_in') && ~isempty(obj.config.disagg_apply_in)
                        use_new_disagg = strcmpi(string(obj.config.disagg_apply_in), "rhs");
                    else
                        use_new_disagg = true;
                    end
                end
            end
            if use_new_disagg && isempty(obj.grid)
                use_new_disagg = false;
            end

            % disagg_rate (per day)
            disagg_rate = 0.0;
            try
                if isprop(obj.config,'disagg_rate') && ~isempty(obj.config.disagg_rate)
                    disagg_rate = obj.config.disagg_rate;
                end
            catch
                disagg_rate = 0.0;
            end
            if ~isfinite(disagg_rate), disagg_rate = 0.0; end
            disagg_rate = max(disagg_rate, 0.0);

            % split state
            v_lin = v;
            if clip_negative
                v_nl = max(v, 0);
            else
                v_nl = v;
            end

            % growth mode
            growth_mode = "shift";
            if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                growth_mode = lower(string(obj.config.growth_mode));
            end
            mu = obj.config.growth;

            % betas
            [b1,b2,b3,b4,b5] = obj.getBetasForRHS();

            % ==========================================================
            % Column branch
            % ==========================================================
            if isprop(obj.config, 'use_column') && obj.config.use_column
                Nz = obj.config.getNumLayers();

                if numel(v) ~= Ns * Nz
                    error('Column RHS: state length mismatch. Expected %d, got %d.', Ns*Nz, numel(v));
                end

                N  = reshape(v_nl, [Ns, Nz]);
                dN = zeros(Ns, Nz);

                % z-grid
                zc = [];
                try
                    if ~isempty(obj.z_cache) && numel(obj.z_cache) == Nz
                        zc = obj.z_cache(:);
                    elseif ismethod(obj.config,'getZ')
                        zc = obj.config.getZ();
                        zc = zc(:);
                        obj.z_cache = zc;
                    end
                catch
                    zc = [];
                end
                if isempty(zc)
                    zc = ((1:Nz) - 0.5) * obj.config.dz;
                    zc = zc(:);
                end

                for k = 1:Nz
                    vk = N(:, k);

                    eps_here = obj.getEpsAtLayer(t, k, zc(k));

                    % new disagg as rate term
                    dv_disagg_k = zeros(size(vk));
                    if use_new_disagg && disagg_rate > 0
                        if ~isempty(eps_here) && isfinite(eps_here) && eps_here > 0

                            eps_ref = eps_here;
                            try
                                if isprop(obj.config,'eps_ref') && ~isempty(obj.config.eps_ref)
                                    eps_ref = obj.config.eps_ref;
                                end
                            catch
                                eps_ref = eps_here;
                            end

                            n_exp = 0.45;
                            if isprop(obj.config,'disagg_n_exp') && ~isempty(obj.config.disagg_n_exp)
                                n_exp = obj.config.disagg_n_exp;
                            end

                            vk_before = max(vk, 0);

                            gridS = obj.getSectionGrid();
                            vk_after  = Disaggregation.applyWithScaling(vk_before, gridS, eps_here, eps_ref, n_exp);

                            dv_disagg_k = disagg_rate * (vk_after - vk_before);
                        end
                    end

                    % coag
                    vk_safe  = max(vk, 0);
                    vk_r     = vk_safe';
                    vk_shift = [0, vk_r(1:Ns-1)];

                    if ~enable_coag
                        term1 = 0 * vk_r;
                        term2 = 0 * vk_r;
                    else
                        term_gain2 = vk_r * b2;  term_gain2 = vk_r .* term_gain2;
                        term_loss3 = vk_r * b3;  term_loss3 = vk_r .* term_loss3;
                        term_loss4 = vk_r * b4;  term_loss4 = vk_r .* term_loss4;
                        term_loss5 = vk_r * b5;  term_loss5 = vk_r .* term_loss5;
                        term1 = (term_gain2 - term_loss3 - term_loss4 - term_loss5);

                        term2 = vk_r * b1;
                        term2 = term2 .* vk_shift;
                    end

                    dvk = (term1 + term2)';

                    % NEW: enforce biovolume closure for coag (per layer)
                    % (REMOVED the old duplicate BV-closure block here)
                    dvk = obj.applyCoagBvClosure(dvk);

                    % legacy matrix disagg (ONLY when enable_disagg==true)
                    if ~use_new_disagg
                        if enable_disagg
                            dvk = dvk - obj.disaggMinus * vk_safe + obj.disaggPlus * vk_safe;
                        end
                    end

                    dvk = dvk + dv_disagg_k;
                    dN(:, k) = dvk;
                end

                % linear operator on RAW state
                term3_flat = obj.linear * v_lin;
                term3 = reshape(term3_flat, [Ns, Nz]);
                dN = dN + term3;

                % legacy growth_mode == "pp"
                if growth_mode == "pp"
                    dN(1,1) = dN(1,1) + mu * v_lin(1);
                end

                % explicit PP source
                dv_pp = obj.getPPSourceVector(v_lin);
                if any(dv_pp ~= 0)
                    dN = dN + reshape(dv_pp, [Ns, Nz]);
                end

                dvdt = dN(:);

                % debug
                try
                    if isprop(obj.config,'debug_rhs_units') && logical(obj.config.debug_rhs_units)
                        dv_lin_dbg = obj.linear * v_lin;
                        dv_rhs_dbg = dvdt;
                        fprintf('[DEBUG RHS UNITS] t=%.4f d | ||A*v||_inf=%.3e | ||rhs||_inf=%.3e | ratio=%.3f\n', ...
                            t, norm(dv_lin_dbg, inf), norm(dv_rhs_dbg, inf), norm(dv_lin_dbg, inf)/max(norm(dv_rhs_dbg, inf), eps));
                    end
                catch
                end

                return;
            end

            % ==========================================================
            % 0-D slab branch
            % ==========================================================
            v_r     = v_nl';
            v_shift = [0, v_r(1:Ns-1)];

            if ~enable_coag
                term1 = 0 * v_r;
                term2 = 0 * v_r;
            else
                term_gain2 = v_r * b2;  term_gain2 = v_r .* term_gain2;
                term_loss3 = v_r * b3;  term_loss3 = v_r .* term_loss3;
                term_loss4 = v_r * b4;  term_loss4 = v_r .* term_loss4;
                term_loss5 = v_r * b5;  term_loss5 = v_r .* term_loss5;
                term1 = (term_gain2 - term_loss3 - term_loss4 - term_loss5);

                term2 = v_r * b1;
                term2 = term2 .* v_shift;
            end

            term3 = obj.linear * v_lin;
            dv_disagg = 0 * v_lin;

            if ~use_new_disagg
                if enable_disagg
                    term4 = -obj.disaggMinus * max(v_nl,0);
                    term5 =  obj.disaggPlus  * max(v_nl,0);
                else
                    term4 = 0 * v_lin;
                    term5 = 0 * v_lin;
                end
            else
                term4 = 0 * v_lin;
                term5 = 0 * v_lin;

                if disagg_rate > 0
                    eps_here = obj.getEpsScalar(t);
                    if ~isempty(eps_here) && isfinite(eps_here) && eps_here > 0
                        eps_ref = eps_here;
                        if isprop(obj.config,'eps_ref') && ~isempty(obj.config.eps_ref)
                            eps_ref = obj.config.eps_ref;
                        end
                        n_exp = 0.45;
                        if isprop(obj.config,'disagg_n_exp') && ~isempty(obj.config.disagg_n_exp)
                            n_exp = obj.config.disagg_n_exp;
                        end
                        v_before = max(v_nl,0);

                        gridS = obj.getSectionGrid();
                        v_after  = Disaggregation.applyWithScaling(v_before, gridS, eps_here, eps_ref, n_exp);

                        dv_disagg = disagg_rate * (v_after - v_before);
                    end
                end
            end

            % NEW: apply closure to the coag part in 0-D
            dv_coag0 = (term1 + term2).';
            dv_coag0 = obj.applyCoagBvClosure(dv_coag0);

            dvdt = dv_coag0 + term3 + term4 + term5 + dv_disagg;

            if growth_mode == "pp"
                dvdt(1) = dvdt(1) + mu * v_lin(1);
            end

            dvdt = dvdt + obj.getPPSourceVector(v_lin);
        end

        % ==========================================================
        % Diagnostics decomposition (UPDATED)
        % ==========================================================
        function varargout = decomposeTerms(obj, t, v)
            v = v(:);

            dv_tot = obj.rhs(t, v);

            v_raw  = v;
            dv_lin = obj.linear * v_raw;

            % PP term (legacy + explicit)
            dv_pp = zeros(size(v));

            growth_mode = "shift";
            try
                if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                    growth_mode = lower(string(obj.config.growth_mode));
                end
            catch
            end

            if growth_mode == "pp"
                mu = obj.config.growth;
                dv_pp(1) = dv_pp(1) + mu * v_raw(1);
            end

            dv_pp = dv_pp + obj.getPPSourceVector(v_raw);

            % switches (same logic as evaluate)
            enable_coag = true;
            try
                if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag)
                    enable_coag = logical(obj.config.enable_coag);
                end
            catch
                enable_coag = true;
            end

            enable_disagg = false;
            try
                if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                    enable_disagg = logical(obj.config.enable_disagg);
                end
            catch
                enable_disagg = false;
            end

            use_new_disagg = false;
            try
                if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg) && logical(obj.config.enable_disagg)
                    if isprop(obj.config,'disagg_apply_in') && ~isempty(obj.config.disagg_apply_in)
                        use_new_disagg = strcmpi(string(obj.config.disagg_apply_in), "rhs");
                    else
                        use_new_disagg = true;
                    end
                end
            catch
                use_new_disagg = false;
            end
            if use_new_disagg && isempty(obj.grid)
                use_new_disagg = false;
            end

            disagg_rate = 0.0;
            try
                if isprop(obj.config,'disagg_rate') && ~isempty(obj.config.disagg_rate)
                    disagg_rate = obj.config.disagg_rate;
                end
            catch
                disagg_rate = 0.0;
            end
            if ~isfinite(disagg_rate), disagg_rate = 0.0; end
            disagg_rate = max(disagg_rate, 0.0);

            clip_negative = true;
            try
                if isprop(obj.config,'clip_negative')
                    clip_negative = logical(obj.config.clip_negative);
                end
            catch
                clip_negative = true;
            end

            if clip_negative
                v_nl = max(v,0);
            else
                v_nl = v;
            end

            % --------------------------
            % Disaggregation term
            % --------------------------
            dv_disagg = zeros(size(v));

            if use_new_disagg && (disagg_rate > 0)
                Ns = obj.config.n_sections;

                if isprop(obj.config,'use_column') && obj.config.use_column
                    Nz = obj.config.getNumLayers();

                    if numel(v) == Ns*Nz
                        N = reshape(max(v,0), [Ns, Nz]);

                        zc = [];
                        try
                            if ~isempty(obj.z_cache) && numel(obj.z_cache) == Nz
                                zc = obj.z_cache(:);
                            elseif ismethod(obj.config,'getZ')
                                zc = obj.config.getZ();
                                zc = zc(:);
                                obj.z_cache = zc;
                            end
                        catch
                            zc = [];
                        end
                        if isempty(zc)
                            zc = ((1:Nz) - 0.5) * obj.config.dz;
                            zc = zc(:);
                        end

                        for k = 1:Nz
                            vk_before = N(:,k);

                            eps_here = obj.getEpsAtLayer(t, k, zc(k));

                            if ~isempty(eps_here) && isfinite(eps_here) && eps_here > 0
                                eps_ref = eps_here;
                                if isprop(obj.config,'eps_ref') && ~isempty(obj.config.eps_ref)
                                    eps_ref = obj.config.eps_ref;
                                end
                                n_exp = 0.45;
                                if isprop(obj.config,'disagg_n_exp') && ~isempty(obj.config.disagg_n_exp)
                                    n_exp = obj.config.disagg_n_exp;
                                end

                                gridS = obj.getSectionGrid();
                                vk_after = Disaggregation.applyWithScaling(vk_before, gridS, eps_here, eps_ref, n_exp);

                                dv_disagg((k-1)*Ns + (1:Ns)) = disagg_rate * (vk_after - vk_before);
                            end
                        end
                    end
                else
                    v_before = max(v,0);
                    eps_here = obj.getEpsScalar(t);

                    if ~isempty(eps_here) && isfinite(eps_here) && eps_here > 0
                        eps_ref = eps_here;
                        if isprop(obj.config,'eps_ref') && ~isempty(obj.config.eps_ref)
                            eps_ref = obj.config.eps_ref;
                        end
                        n_exp = 0.45;
                        if isprop(obj.config,'disagg_n_exp') && ~isempty(obj.config.disagg_n_exp)
                            n_exp = obj.config.disagg_n_exp;
                        end

                        gridS = obj.getSectionGrid();
                        v_after  = Disaggregation.applyWithScaling(v_before, gridS, eps_here, eps_ref, n_exp);

                        dv_disagg = disagg_rate * (v_after - v_before);
                    end
                end

            elseif ~use_new_disagg
                if enable_disagg
                    v_nonneg = max(v,0);
                    dv_disagg = (-obj.disaggMinus * v_nonneg) + (obj.disaggPlus * v_nonneg);
                else
                    dv_disagg = zeros(size(v));
                end
            end

            % --------------------------
            % Coagulation term
            % --------------------------
            dv_coag = zeros(size(v));

            if enable_coag
                Ns = obj.config.n_sections;
                [b1,b2,b3,b4,b5] = obj.getBetasForRHS();

                if isprop(obj.config,'use_column') && obj.config.use_column
                    Nz = obj.config.getNumLayers();

                    if numel(v) == Ns*Nz
                        N = reshape(v_nl, [Ns, Nz]);
                        dNcoag = zeros(Ns, Nz);

                        for k = 1:Nz
                            vk = N(:,k);
                            vk_safe  = max(vk, 0);
                            vk_r     = vk_safe';
                            vk_shift = [0, vk_r(1:Ns-1)];

                            term_gain2 = vk_r * b2;  term_gain2 = vk_r .* term_gain2;
                            term_loss3 = vk_r * b3;  term_loss3 = vk_r .* term_loss3;
                            term_loss4 = vk_r * b4;  term_loss4 = vk_r .* term_loss4;
                            term_loss5 = vk_r * b5;  term_loss5 = vk_r .* term_loss5;
                            term1 = (term_gain2 - term_loss3 - term_loss4 - term_loss5);

                            term2 = vk_r * b1;
                            term2 = term2 .* vk_shift;

                            tmp = (term1 + term2).';
                            tmp = obj.applyCoagBvClosure(tmp);
                            dNcoag(:,k) = tmp;
                        end

                        dv_coag = dNcoag(:);
                    end
                else
                    v_r     = v_nl';
                    v_shift = [0, v_r(1:Ns-1)];

                    term_gain2 = v_r * b2;  term_gain2 = v_r .* term_gain2;
                    term_loss3 = v_r * b3;  term_loss3 = v_r .* term_loss3;
                    term_loss4 = v_r * b4;  term_loss4 = v_r .* term_loss4;
                    term_loss5 = v_r * b5;  term_loss5 = v_r .* term_loss5;
                    term1 = (term_gain2 - term_loss3 - term_loss4 - term_loss5);

                    term2 = v_r * b1;
                    term2 = term2 .* v_shift;

                    dv_coag = (term1 + term2).';
                    dv_coag = obj.applyCoagBvClosure(dv_coag);
                end
            end

            dv_other = dv_tot - (dv_lin + dv_pp + dv_coag + dv_disagg);

            S = struct();
            S.dv_tot    = dv_tot;
            S.dv_lin    = dv_lin;
            S.dv_pp     = dv_pp;
            S.dv_coag   = dv_coag;
            S.dv_disagg = dv_disagg;
            S.dv_other  = dv_other;

            if nargout <= 1
                varargout{1} = S;
                return;
            end

            varargout{1} = dv_tot;
            if nargout >= 2, varargout{2} = dv_lin;    end
            if nargout >= 3, varargout{3} = dv_coag;   end
            if nargout >= 4, varargout{4} = dv_pp;     end
            if nargout >= 5, varargout{5} = dv_disagg; end
            if nargout >= 6, varargout{6} = dv_other;  end
        end

        function J = jacobian(obj, t, v) %#ok<INUSD>
            use_full_jacobian = false;
            if isprop(obj.config,'use_full_jacobian') && ~isempty(obj.config.use_full_jacobian)
                use_full_jacobian = logical(obj.config.use_full_jacobian);
            end

            enable_disagg = false;
            try
                if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                    enable_disagg = logical(obj.config.enable_disagg);
                end
            catch
                enable_disagg = false;
            end

            if isprop(obj.config, 'use_column') && obj.config.use_column
                Ns = obj.config.n_sections;
                Nz = obj.config.getNumLayers();
                Iz = speye(Nz);

                use_new_disagg = false;
                if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                    if logical(obj.config.enable_disagg)
                        if isprop(obj.config,'disagg_apply_in') && ~isempty(obj.config.disagg_apply_in)
                            use_new_disagg = strcmpi(string(obj.config.disagg_apply_in), "rhs");
                        else
                            use_new_disagg = true;
                        end
                    end
                end
                if use_new_disagg && isempty(obj.grid)
                    use_new_disagg = false;
                end

                if use_new_disagg
                    D = sparse(Ns, Ns);
                else
                    if enable_disagg
                        D = (-obj.disaggMinus + obj.disaggPlus);
                    else
                        D = sparse(Ns, Ns);
                    end
                end

                J = obj.linear + kron(Iz, D);

                growth_mode = "shift";
                if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                    growth_mode = lower(string(obj.config.growth_mode));
                end
                if growth_mode == "pp"
                    mu = obj.config.growth;
                    J(1,1) = J(1,1) + mu;
                end

                % explicit PP source is additive constant => zero Jacobian
                return;
            end

            if ~use_full_jacobian
                if enable_disagg
                    J = obj.linear - obj.disaggMinus + obj.disaggPlus;
                else
                    J = obj.linear;
                end
                return;
            end

            % Full Jacobian path (kept)
            try
                [b1,b2,b3,b4,b5] = obj.getBetasForRHS(); %#ok<ASGLU>
                n_sections = length(v);
                v_r = v';

                v_mat   = v_r(ones(1, n_sections), :);
                v_shift = [zeros(n_sections, 1), v_mat(:, 1:end-1)];

                bnet = (b2 - b3 - b4 - b5);

                term1_tmp = v_r * bnet;
                term1     = diag(term1_tmp) + diag(v_r) * bnet;

                term2a = v_r * b1;
                term2a = diag(term2a(2:end), -1);

                term2b = diag(b1, 1);
                term2b = term2b' .* v_r(1:end-1);
                term2b = diag(term2b, -1);

                term2c = diag(v_r(2:end), -1) .* bnet';
                term2  = term2a + term2b + term2c;

                term3a = b1  .* v_shift;
                term3b = bnet .* v_mat;

                term3  = (term3a + term3b)';
                term3  = triu(term3, 2) + tril(term3, -1);

                if enable_disagg
                    J = term1 + term2 + term3 + obj.linear - obj.disaggMinus + obj.disaggPlus;
                else
                    J = term1 + term2 + term3 + obj.linear;
                end

            catch ME
                warning('Full Jacobian failed (%s). Falling back to safe Jacobian.', ME.message);

                if enable_disagg
                    J = obj.linear - obj.disaggMinus + obj.disaggPlus;
                else
                    J = obj.linear;
                end

                growth_mode = "shift";
                if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                    growth_mode = lower(string(obj.config.growth_mode));
                end
                if growth_mode == "pp"
                    mu = obj.config.growth;
                    J(1,1) = J(1,1) + mu;
                end
                return;
            end
        end

        function [term1, term2, term3, term4, term5] = rateTerms(obj, v)
            % kept as your existing routine (linear + legacy disagg matrices)
            % NOTE: Part-2 uses decomposeTerms(), not rateTerms().

            Ns = obj.config.n_sections;

            enable_coag = true;
            if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag)
                enable_coag = logical(obj.config.enable_coag);
            end

            enable_disagg = false;
            try
                if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                    enable_disagg = logical(obj.config.enable_disagg);
                end
            catch
                enable_disagg = false;
            end

            clip_negative = true;
            if isprop(obj.config, 'clip_negative')
                clip_negative = logical(obj.config.clip_negative);
            end

            use_new_disagg = false;
            if isprop(obj.config,'enable_disagg') && ~isempty(obj.config.enable_disagg)
                if logical(obj.config.enable_disagg)
                    if isprop(obj.config,'disagg_apply_in') && ~isempty(obj.config.disagg_apply_in)
                        use_new_disagg = strcmpi(string(obj.config.disagg_apply_in), "rhs");
                    else
                        use_new_disagg = true;
                    end
                end
            end
            if use_new_disagg && isempty(obj.grid)
                use_new_disagg = false;
            end

            v_lin = v;
            if clip_negative
                v_nl = max(v,0);
            else
                v_nl = v;
            end

            has_b25 = false;
            try
                has_b25 = isprop(obj.betas,'b25') && ~isempty(obj.betas.b25);
            catch
                has_b25 = false;
            end

            [b1,b2,b3,b4,b5] = obj.getBetasForRHS();

            if isprop(obj.config, 'use_column') && obj.config.use_column
                Nz = obj.config.getNumLayers();

                if numel(v) ~= Ns * Nz
                    error('rateTerms: Column state length mismatch. Expected %d, got %d.', Ns*Nz, numel(v));
                end

                N = reshape(v_nl, [Ns, Nz]);

                T1 = zeros(Ns, Nz);
                T2 = zeros(Ns, Nz);

                for k = 1:Nz
                    vk = N(:, k);
                    vk_safe = max(vk, 0);

                    vk_r     = vk_safe';
                    vk_shift = [0, vk_r(1:Ns-1)];

                    if ~enable_coag
                        a1 = 0 * vk_r;
                        a2 = 0 * vk_r;
                    else
                        if has_b25
                            a1 = vk_r * obj.betas.b25;
                            a1 = vk_r .* a1;
                        else
                            bnet = (b2 - b3 - b4 - b5);
                            a1 = vk_r * bnet;
                            a1 = vk_r .* a1;
                        end

                        a2 = vk_r * b1;
                        a2 = a2 .* vk_shift;
                    end

                    T1(:, k) = a1';
                    T2(:, k) = a2';
                end

                term1 = T1(:);
                term2 = T2(:);

                term3 = obj.linear * v_lin;

                Dm = obj.disaggMinus;
                Dp = obj.disaggPlus;

                T4 = zeros(Ns, Nz);
                T5 = zeros(Ns, Nz);

                for k = 1:Nz
                    vk = N(:, k);
                    vk_safe = max(vk, 0);

                    if ~use_new_disagg
                        if enable_disagg
                            T4(:, k) = -Dm * vk_safe;
                            T5(:, k) =  Dp * vk_safe;
                        else
                            T4(:, k) = 0 * vk_safe;
                            T5(:, k) = 0 * vk_safe;
                        end
                    else
                        T4(:, k) = 0 * vk_safe;
                        T5(:, k) = 0 * vk_safe;
                    end
                end

                term4 = T4(:);
                term5 = T5(:);
                return;
            end

            % 0-D
            v_r     = v_nl';
            v_shift = [0, v_r(1:Ns-1)];

            if ~enable_coag
                term1 = 0 * v_r;
                term2 = 0 * v_r;
            else
                if has_b25
                    term1 = v_r * obj.betas.b25;
                    term1 = v_r .* term1;
                else
                    bnet = (b2 - b3 - b4 - b5);
                    term1 = v_r * bnet;
                    term1 = v_r .* term1;
                end

                term2 = v_r * b1;
                term2 = term2 .* v_shift;
            end

            term3 = obj.linear * v_lin;

            if ~use_new_disagg
                if enable_disagg
                    term4 = -obj.disaggMinus * max(v_nl,0);
                    term5 =  obj.disaggPlus  * max(v_nl,0);
                else
                    term4 = 0 * v_lin;
                    term5 = 0 * v_lin;
                end
            else
                term4 = 0 * v_lin;
                term5 = 0 * v_lin;
            end
        end

        function validate(obj)
            if isempty(obj.betas) || isempty(obj.linear)
                error('RHS not properly initialized');
            end

            Ns = obj.config.n_sections;

            if isprop(obj.config, 'use_column') && obj.config.use_column
                Nz = obj.config.getNumLayers();
                expected = Ns * Nz;
                if size(obj.linear, 1) ~= expected
                    error('Linear matrix dimension mismatch (column). Expected %d, got %d.', expected, size(obj.linear, 1));
                end
            else
                if size(obj.linear, 1) ~= Ns
                    error('Linear matrix dimension mismatch (0-D). Expected %d, got %d.', Ns, size(obj.linear, 1));
                end
            end

            if obj.betas.getNumSections() ~= Ns
                error('Beta matrices dimension mismatch');
            end
        end
    end
end
