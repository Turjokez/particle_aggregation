% ==============================================================
% FILE: CoagulationRHS.m
% ==============================================================
% SAFE VERSION:
% - b25 is treated as DIAGNOSTIC by default (NOT used in RHS)
% - legacy b25 path is kept but guarded by cfg.use_legacy_b25=true
%   AND a sanity check (finite + nonnegative)
% ==============================================================

classdef CoagulationRHS < handle

    properties
        betas
        linear
        disaggMinus
        disaggPlus
        config

        eps_time  = [];
        eps_vals  = [];
        eps_const = [];

        grid      = [];
        z_cache   = [];

        printed_b25_stats = false; % NEW
    end

    methods
        function obj = CoagulationRHS(betas, linear, disaggMinus, disaggPlus, config, varargin)
            obj.betas       = betas;
            obj.linear      = linear;
            obj.disaggMinus = disaggMinus;
            obj.disaggPlus  = disaggPlus;
            obj.config      = config;

            if isprop(config,'epsilon_const')
                obj.eps_const = config.epsilon_const;
            else
                obj.eps_const = [];
            end

            if nargin >= 6 && ~isempty(varargin)
                try, obj.grid = varargin{1}; catch, obj.grid = []; end
            end

            try
                if ismethod(config,'getZ')
                    obj.z_cache = config.getZ();
                end
            catch
                obj.z_cache = [];
            end
        end

        function setEpsilonSeries(obj, t, eps_vals)
            obj.eps_time = t(:);

            if isvector(eps_vals)
                obj.eps_vals = eps_vals(:);
                n = min(numel(obj.eps_time), numel(obj.eps_vals));
                obj.eps_time = obj.eps_time(1:n);
                obj.eps_vals = obj.eps_vals(1:n);
            else
                try
                    Nt = numel(obj.eps_time);
                    if size(eps_vals,2) == Nt && size(eps_vals,1) ~= Nt
                        eps_vals = eps_vals.';
                    end
                catch
                end
                obj.eps_vals = eps_vals;
            end
        end

        function setGrid(obj, gridObj)
            obj.grid = gridObj;
        end

        function dvdt = rhs(obj, t, v)
            dvdt = obj.evaluate(t, v);
        end

        function dvdt = evaluate(obj, t, v)
            Ns = obj.config.n_sections;

            enable_coag = true;
            if isprop(obj.config,'enable_coag') && ~isempty(obj.config.enable_coag)
                enable_coag = logical(obj.config.enable_coag);
            end

            clip_negative = true;
            if isprop(obj.config,'clip_negative') && ~isempty(obj.config.clip_negative)
                clip_negative = logical(obj.config.clip_negative);
            end

            [use_new_disagg, disagg_rate] = obj.getDisaggControls();
            [enable_pp, pp_rate, pp_bin, pp_layer] = obj.getPPControls(Ns);

            enforce_closure = false;
            if isprop(obj.config,'enforce_coag_bv_closure') && ~isempty(obj.config.enforce_coag_bv_closure)
                enforce_closure = logical(obj.config.enforce_coag_bv_closure);
            end

            growth_mode = "shift";
            if isprop(obj.config,'growth_mode') && ~isempty(obj.config.growth_mode)
                growth_mode = lower(string(obj.config.growth_mode));
            end
            mu = obj.config.growth;

            state_is_bv = false;
            if isprop(obj.config,'state_is_biovolume') && ~isempty(obj.config.state_is_biovolume)
                state_is_bv = logical(obj.config.state_is_biovolume);
            end

            [av_vol, have_av] = obj.getAvVolForCoag(Ns);
            if state_is_bv && ~have_av
                error('CoagulationRHS:MissingAvVol', 'state_is_biovolume=true requires grid.av_vol in RHS.');
            end

            [b1, b2, b3, b4, b5, b25] = obj.getBetasRowForm();

            % ==========================================================
            % NEW-2026-01-14: b25 is DIAGNOSTIC by default
            % Guard legacy b25 usage behind cfg.use_legacy_b25
            % ==========================================================
            use_legacy_b25 = false;
            try
                if isprop(obj.config,'use_legacy_b25') && ~isempty(obj.config.use_legacy_b25)
                    use_legacy_b25 = logical(obj.config.use_legacy_b25);
                end
            catch
                use_legacy_b25 = false;
            end

            has_b25 = ~isempty(b25);
            b25_ok = false;
            if has_b25
                try
                    b25_ok = all(isfinite(b25(:))) && (min(b25(:)) >= 0);
                catch
                    b25_ok = false;
                end
            end

            use_b25_path = enable_coag && use_legacy_b25 && has_b25 && b25_ok;

            % debug print at t0
            try
                do_dbg = false;
                if isprop(obj.config,'debug_coag_leak') && ~isempty(obj.config.debug_coag_leak)
                    do_dbg = logical(obj.config.debug_coag_leak);
                end
                if do_dbg && ~obj.printed_b25_stats && (abs(t) < 1e-12 || abs(t-obj.config.t_init) < 1e-12)
                    if has_b25
                        fprintf('[b25 stats] min=%.3e max=%.3e | b25_ok=%d | use_legacy_b25=%d\n', ...
                            min(b25(:)), max(b25(:)), b25_ok, use_legacy_b25);
                    else
                        fprintf('[b25 stats] b25 not present\n');
                    end
                    obj.printed_b25_stats = true;
                end
            catch
            end

            % ==========================================================
            % COLUMN MODE
            % ==========================================================
            if isprop(obj.config,'use_column') && ~isempty(obj.config.use_column) && obj.config.use_column
                Nz = obj.config.getNumLayers();

                if numel(v) ~= Ns*Nz
                    error('Column RHS: state length mismatch. Expected %d, got %d.', Ns*Nz, numel(v));
                end

                v_lin = v(:);
                v_nl  = v_lin;
                if clip_negative, v_nl = max(v_nl,0); end

                V  = reshape(v_nl, [Ns, Nz]);
                dV = zeros(Ns, Nz);

                zc = obj.getZCenters(Nz);

                has_eps_fun = false;
                try
                    has_eps_fun = isprop(obj.config,'eps_fun') && ~isempty(obj.config.eps_fun) && isa(obj.config.eps_fun,'function_handle');
                catch
                    has_eps_fun = false;
                end

                has_eps_matrix = false;
                try
                    if isprop(obj.config,'epsilon_time') && isprop(obj.config,'epsilon_series')
                        if ~isempty(obj.config.epsilon_time) && ~isempty(obj.config.epsilon_series)
                            if ~isvector(obj.config.epsilon_series) && size(obj.config.epsilon_series,2) == Nz
                                has_eps_matrix = true;
                            end
                        end
                    end
                catch
                    has_eps_matrix = false;
                end

                for k = 1:Nz
                    vk = max(V(:,k),0);

                    if state_is_bv
                        nk = vk ./ max(av_vol, eps);
                    else
                        nk = vk;
                    end

                    nr = nk.';
                    nsh = [0, nr(1:Ns-1)];

                    if ~enable_coag
                        dn = zeros(Ns,1);
                    else
                        if use_b25_path
                            % ---- LEGACY b25 path (kept but guarded) ----
                            term1 = nr * b25; term1 = nr .* term1;
                            term2 = nr * b1;  term2 = term2 .* nsh;
                            dn = (term1 + term2).';
                        else
                            % ---- SAFE default path ----
                            g2 = nr * b2; g2 = nr .* g2;
                            l3 = nr * b3; l3 = nr .* l3;
                            l4 = nr * b4; l4 = nr .* l4;
                            l5 = nr * b5; l5 = nr .* l5;

                            term1 = (g2 - l3 - l4 - l5);
                            term2 = nr * b1; term2 = term2 .* nsh;
                            dn = (term1 + term2).';
                        end
                    end

                    obj.debugCoagLeak(t, dn, av_vol, sprintf('COL k=%d',k));

                    if enforce_closure && enable_coag
                        leakM = sum(dn .* av_vol);
                        dn(end) = dn(end) - leakM / av_vol(end);
                    end

                    if state_is_bv
                        dvk = dn .* av_vol;
                    else
                        dvk = dn;
                    end

                    if ~use_new_disagg
                        dvk = dvk - obj.disaggMinus*vk + obj.disaggPlus*vk;
                    elseif disagg_rate > 0
                        eps_here = obj.getEpsHereColumn(t, k, zc, Nz, has_eps_fun, has_eps_matrix);
                        dvk = dvk + obj.computeNewDisaggTendencyVec(t, vk, eps_here, disagg_rate, k);
                    end

                    dV(:,k) = dvk;
                end

                % linear (BIOVOLUME-based)
                if state_is_bv
                    dv_lin_flat = obj.linear * v_lin;
                else
                    av_rep = repmat(av_vol(:), Nz, 1);
                    v_bv   = v_lin .* av_rep;
                    dv_bv  = obj.linear * v_bv;
                    dv_lin_flat = dv_bv ./ max(av_rep, eps);
                end

                dV = dV + reshape(dv_lin_flat, [Ns, Nz]);

                if growth_mode == "pp"
                    dV(1,1) = dV(1,1) + mu * v_lin(1);
                end

                if enable_pp && pp_rate ~= 0
                    kb = max(1, min(Ns, pp_bin));
                    kz = max(1, min(Nz, pp_layer));
                    dV(kb,kz) = dV(kb,kz) + pp_rate;
                end

                dvdt = dV(:);
                return;
            end

            % ==========================================================
            % 0-D MODE
            % ==========================================================
            v_lin = v(:);
            v_nl  = v_lin;
            if clip_negative, v_nl = max(v_nl,0); end

            if state_is_bv
                n = v_nl ./ max(av_vol, eps);
            else
                n = v_nl;
            end

            nr  = n.';
            nsh = [0, nr(1:Ns-1)];

            if ~enable_coag
                dn = zeros(Ns,1);
            else
                if use_b25_path
                    term1 = nr * b25; term1 = nr .* term1;
                    term2 = nr * b1;  term2 = term2 .* nsh;
                    dn = (term1 + term2).';
                else
                    g2 = nr * b2; g2 = nr .* g2;
                    l3 = nr * b3; l3 = nr .* l3;
                    l4 = nr * b4; l4 = nr .* l4;
                    l5 = nr * b5; l5 = nr .* l5;

                    term1 = (g2 - l3 - l4 - l5);
                    term2 = nr * b1; term2 = term2 .* nsh;
                    dn = (term1 + term2).';
                end
            end

            obj.debugCoagLeak(t, dn, av_vol, '0D');

            if enforce_closure && enable_coag
                leakM = sum(dn .* av_vol);
                dn(end) = dn(end) - leakM / av_vol(end);
            end

            if state_is_bv
                dv_coag = dn .* av_vol;
                dv_lin  = obj.linear * v_lin;
            else
                dv_coag = dn;
                v_bv    = v_lin .* av_vol;
                dv_bv   = obj.linear * v_bv;
                dv_lin  = dv_bv ./ max(av_vol, eps);
            end

            if ~use_new_disagg
                term4 = -obj.disaggMinus * v_nl;
                term5 =  obj.disaggPlus  * v_nl;
                dv_disagg = 0*v_lin;
            else
                term4 = 0*v_lin; term5 = 0*v_lin;
                dv_disagg = 0*v_lin;
                if disagg_rate > 0
                    eps_here = obj.getEpsScalar(t);
                    dv_disagg = obj.computeNewDisaggTendencyVec(t, v_nl, eps_here, disagg_rate, []);
                end
            end

            dvdt = dv_coag + dv_lin + term4 + term5 + dv_disagg;

            if growth_mode == "pp"
                dvdt(1) = dvdt(1) + mu * v_lin(1);
            end

            if enable_pp && pp_rate ~= 0
                ib = max(1, min(Ns, pp_bin));
                dvdt(ib) = dvdt(ib) + pp_rate;
            end
        end
    end

    % ==========================================================
    % Helpers (unchanged structure)
    % ==========================================================
    methods (Access=private)

        function debugCoagLeak(obj, t, dn, av_vol, tag)
            try
                do_dbg = false;
                if isprop(obj.config,'debug_coag_leak') && ~isempty(obj.config.debug_coag_leak)
                    do_dbg = logical(obj.config.debug_coag_leak);
                end
                if do_dbg
                    leakM = sum(dn .* av_vol);
                    if abs(t) < 1e-12 || abs(t-obj.config.t_init) < 1e-12
                        fprintf('[COAG DEBUG] %s t=%.6g leak(sum(dn.*av)) = %.6e\n', tag, t, leakM);
                    end
                end
            catch
            end
        end

        function [b1,b2,b3,b4,b5,b25] = getBetasRowForm(obj)
            b1 = obj.betas.b1;
            b2 = obj.betas.b2;
            b3 = obj.betas.b3;
            b4 = obj.betas.b4;
            b5 = obj.betas.b5;

            b25 = [];
            try
                if isprop(obj.betas,'b25') && ~isempty(obj.betas.b25)
                    b25 = obj.betas.b25;
                end
            catch
                b25 = [];
            end

            do_fix = false;
            if isprop(obj.config,'beta_row_fix') && ~isempty(obj.config.beta_row_fix)
                do_fix = logical(obj.config.beta_row_fix);
            end
            if ~do_fix, return; end

            do_t = false;
            if isprop(obj.config,'beta_transpose_b1b3') && ~isempty(obj.config.beta_transpose_b1b3)
                do_t = logical(obj.config.beta_transpose_b1b3);
            end

            if do_t
                b1 = b1.';
                b3 = b3.';
                if ~isempty(b25), b25 = b25.'; end
            end
        end

        function [use_new_disagg, disagg_rate] = getDisaggControls(obj)
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
        end

        function [enable_pp, pp_rate, pp_bin, pp_layer] = getPPControls(obj, Ns)
            enable_pp = false; pp_rate = 0.0; pp_bin = 1; pp_layer = 1;

            try
                if isprop(obj.config,'enable_pp') && ~isempty(obj.config.enable_pp)
                    enable_pp = logical(obj.config.enable_pp);
                end
            catch
                enable_pp = false;
            end

            try
                if isprop(obj.config,'pp_rate') && ~isempty(obj.config.pp_rate)
                    pp_rate = obj.config.pp_rate;
                end
            catch
                pp_rate = 0.0;
            end
            if ~isfinite(pp_rate), pp_rate = 0.0; end

            try
                if isprop(obj.config,'pp_bin') && ~isempty(obj.config.pp_bin)
                    pp_bin = obj.config.pp_bin;
                end
            catch
                pp_bin = 1;
            end
            pp_bin = max(1, min(Ns, round(pp_bin)));

            try
                if isprop(obj.config,'pp_layer') && ~isempty(obj.config.pp_layer)
                    pp_layer = obj.config.pp_layer;
                end
            catch
                pp_layer = 1;
            end
            pp_layer = max(1, round(pp_layer));
        end

        function eps_now = getEpsScalar(obj, t)
            eps_now = [];
            if ~isempty(obj.eps_time) && ~isempty(obj.eps_vals) && isvector(obj.eps_vals)
                n  = min(numel(obj.eps_time), numel(obj.eps_vals));
                tt = obj.eps_time(1:n);
                ee = obj.eps_vals(1:n);
                try
                    eps_now = interp1(tt(:), ee(:), t, 'linear', 'extrap');
                catch
                    eps_now = [];
                end
            elseif ~isempty(obj.eps_const)
                eps_now = obj.eps_const;
            end
            try
                if ~isempty(eps_now) && isfinite(eps_now)
                    eps_now = max(eps_now, 0);
                end
            catch
            end
        end

        function zc = getZCenters(obj, Nz)
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
        end

        function eps_here = getEpsHereColumn(obj, t, k, zc, Nz, has_eps_fun, has_eps_matrix)
            eps_here = [];
            if has_eps_fun
                try, eps_here = obj.config.eps_fun(t, zc(k)); catch, eps_here = []; end
            end
            if isempty(eps_here) && has_eps_matrix
                try
                    tt = obj.config.epsilon_time(:);
                    Em = obj.config.epsilon_series;
                    if size(Em,2) == Nz
                        eps_here = interp1(tt, Em(:,k), t, 'linear', 'extrap');
                    end
                catch
                    eps_here = [];
                end
            end
            if isempty(eps_here)
                eps_here = obj.getEpsScalar(t);
            end
        end

        function [av, have_av] = getAvVolForCoag(obj, Ns)
            av = []; have_av = false;
            try
                if ~isempty(obj.grid) && isprop(obj.grid,'av_vol') && ~isempty(obj.grid.av_vol)
                    av = obj.grid.av_vol(:);
                end
            catch
                av = [];
            end
            if isempty(av), return; end
            if numel(av) ~= Ns, av = []; return; end
            if any(~isfinite(av)) || any(av <= 0), av = []; return; end
            have_av = true;
        end

        function gridS = getSectionGrid(obj)
            gridS = obj.grid;
            if isempty(gridS)
                error('CoagulationRHS:getSectionGrid:EmptyGrid', ...
                    'obj.grid is empty. Pass grid into CoagulationRHS or call rhs.setGrid(grid).');
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
                error('CoagulationRHS:getSectionGrid:InvalidGrid', 'No radii fallback possible.');
            end
            gridS = SectionGridAdapter(rv, obj.config);
        end

        function dv_disagg = computeNewDisaggTendencyVec(obj, t, vk_before, eps_here, disagg_rate, k) %#ok<INUSD>
            dv_disagg = zeros(size(vk_before));
            if disagg_rate <= 0, return; end
            if isempty(eps_here) || ~isfinite(eps_here) || eps_here <= 0, return; end

            eps_ref = 1e-6;
            try
                if isprop(obj.config,'eps_ref') && ~isempty(obj.config.eps_ref) && obj.config.eps_ref > 0
                    eps_ref = obj.config.eps_ref;
                end
            catch
            end

            n_exp = 0.45;
            if isprop(obj.config,'disagg_n_exp') && ~isempty(obj.config.disagg_n_exp)
                n_exp = obj.config.disagg_n_exp;
            end

            vk_before = max(vk_before, 0);

            gridS = obj.getSectionGrid();
            vk_after = Disaggregation.applyWithScaling(vk_before, gridS, eps_here, eps_ref, n_exp);
            dv_disagg = disagg_rate * (vk_after - vk_before);
        end
    end
end
