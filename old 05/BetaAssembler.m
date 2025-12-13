classdef BetaAssembler < handle
    % =====================================================================
    %  BETAASSEMBLER
    %  Computes sectionally integrated coagulation kernel matrices (b1..b5)
    %  from the selected kernel in SimulationConfig, using the derived grid.
    %
    %  IMPORTANT CONVENTION (must match CoagulationRHS + BetaMatrices):
    %      b25 = b2 + b5   (total loss from section i)
    %
    % =====================================================================

    properties
        config;         % SimulationConfig object (kernel choice, constants)
        grid;           % DerivedGrid object (section limits/coefficients)
        kernels;        % KernelLibrary reference (kept for parity)
    end

    % =====================================================================
    %  PUBLIC API
    % =====================================================================
    methods
        % -----------------------------------------------------------------
        function obj = BetaAssembler(config, grid)
            % Constructor: store references
            obj.config  = config;
            obj.grid    = grid;
            obj.kernels = KernelLibrary(); %#ok<NASGU> % loaded for parity
        end

        % -----------------------------------------------------------------
        function betas = computeFor(obj, kernelName)
            % computeFor: compute matrices for a specific kernel by name.
            % Temporarily overrides config.kernel, computes, then restores.
            original_kernel = obj.config.kernel;
            obj.config.kernel = kernelName;
            betas = obj.computeBetaMatrices();
            obj.config.kernel = original_kernel;
        end

        % -----------------------------------------------------------------
        function betas = combineAndScale(obj, b_brown, b_shear, b_ds)
            % combineAndScale
            %  Combine contributions from Brownian, shear, and differential
            %  sedimentation into a single BetaMatrices, applying consistent
            %  per-day units.

            % ---- Scale each block if present ----------------------------
            if ~isempty(b_brown) && ~isempty(b_brown.b1)
                b_brown_s = obj.scaleBetas(b_brown, obj.grid.conBr * obj.config.day_to_sec);
            else
                b_brown_s = BetaMatrices();
            end

            if ~isempty(b_shear) && ~isempty(b_shear.b1)
                b_shear_s = obj.scaleBetas(b_shear, obj.config.gamma * obj.config.day_to_sec);
            else
                b_shear_s = BetaMatrices();
            end

            if ~isempty(b_ds) && ~isempty(b_ds.b1)
                b_ds_s = obj.scaleBetas(b_ds, obj.grid.setcon * obj.config.day_to_sec);
            else
                b_ds_s = BetaMatrices();
            end

            % ---- Determine target matrix size N from first non-empty ----
            N = 0;
            if ~isempty(b_brown_s.b1), N = size(b_brown_s.b1,1); end
            if N==0 && ~isempty(b_shear_s.b1), N = size(b_shear_s.b1,1); end
            if N==0 && ~isempty(b_ds_s.b1),    N = size(b_ds_s.b1,1);    end

            % ---- Assemble combined matrices (null-protected) -----------
            betas = BetaMatrices();
            if N > 0
                betas.b1 = nz(b_brown_s.b1) + nz(b_shear_s.b1) + nz(b_ds_s.b1);
                betas.b2 = nz(b_brown_s.b2) + nz(b_shear_s.b2) + nz(b_ds_s.b2);
                betas.b3 = nz(b_brown_s.b3) + nz(b_shear_s.b3) + nz(b_ds_s.b3);
                betas.b4 = nz(b_brown_s.b4) + nz(b_shear_s.b4) + nz(b_ds_s.b4);
                betas.b5 = nz(b_brown_s.b5) + nz(b_shear_s.b5) + nz(b_ds_s.b5);

                % >>> CONSISTENT RHS convention: total loss from i
                betas.b25 = betas.b2 + betas.b5;
            end

            % ---- Nested helper: empty → zeros(N) -----------------------
            function A = nz(M)
                if isempty(M), A = zeros(N); else, A = M; end
            end
        end

        % -----------------------------------------------------------------
        function betas_scaled = scaleBetas(~, betas, scale_factor)
            % scaleBetas: multiply each matrix by scale_factor; recompute b25
            betas_scaled      = BetaMatrices();
            betas_scaled.b1   = betas.b1  * scale_factor;
            betas_scaled.b2   = betas.b2  * scale_factor;
            betas_scaled.b3   = betas.b3  * scale_factor;
            betas_scaled.b4   = betas.b4  * scale_factor;
            betas_scaled.b5   = betas.b5  * scale_factor;

            % keep b25 consistent: total loss from i
            betas_scaled.b25  = betas_scaled.b2 + betas_scaled.b5;
        end
    end

    % =====================================================================
    %  PRIVATE IMPLEMENTATION (integration logic unchanged)
    % =====================================================================
    methods (Access = private)

        function betas = computeBetaMatrices(obj)
            % computeBetaMatrices
            %  Port of CalcBetas-style construction of b1..b5.

            n_sections = obj.config.n_sections;
            mlo        = obj.grid.v_lower;  % lower bound of each mass/volume section

            % Initialize all beta matrices with zeros
            beta_init = zeros(n_sections, n_sections);

            % Parameters passed into kernel wrappers
            int_param.amfrac    = obj.grid.amfrac;       % fractal prefactor
            int_param.bmfrac    = obj.grid.bmfrac;       % fractal exponent
            int_param.kernel    = obj.config.kernel;     % selected kernel name
            int_param.r_to_rg   = obj.config.r_to_rg;    % interaction→gyration
            int_param.setcon    = obj.grid.setcon;       % settling constant
            int_param.constants = obj.config;            % pass-through constants

            % -------------------- Case 5: b5 (loss: i & j > i) ----------
            b5 = beta_init;
            for jcol = 1:(n_sections - 1)
                for irow = (jcol + 1):n_sections
                    mj_lo = mlo(jcol);  mj_up = 2.0 * mj_lo;
                    mi_lo = mlo(irow);  mi_up = 2.0 * mi_lo;

                    bndry.mi_lo = mi_lo;   bndry.mi_up = mi_up;
                    bndry.mj_lo = mj_lo;   bndry.mj_up = mj_up;
                    bndry.mjj   = [];      bndry.rjj   = [];
                    bndry.rvjj  = [];

                    b5(irow, jcol) = quadl(@(x) obj.integr5a(x, int_param, bndry), ...
                                           mi_lo, mi_up) / (mi_lo * mj_lo);
                end
            end

            % -------------------- Case 4: b4 (loss: j with itself) -----
            b4 = beta_init;
            for jcol = 1:n_sections
                mj_lo = mlo(jcol);  mj_up = 2.0 * mj_lo;
                mi_lo = mlo(jcol);  mi_up = 2.0 * mi_lo;

                bndry.mi_lo = mi_lo;   bndry.mi_up = mi_up;
                bndry.mj_lo = mj_lo;   bndry.mj_up = mj_up;
                bndry.mjj   = [];      bndry.rjj   = [];
                bndry.rvjj  = [];

                b4(jcol, jcol) = quadl(@(x) obj.integr4a(x, int_param, bndry), ...
                                       mi_lo, mi_up) / (mi_lo * mj_lo);
            end
            b4 = b4 / 2;   % identical-pair double-count correction

            % -------------------- Case 3: b3 (loss: i < j) --------------
            b3 = beta_init;
            for jcol = 2:n_sections
                for irow = 1:(jcol - 1)
                    mj_lo = mlo(jcol);  mj_up = 2.0 * mj_lo;
                    mi_lo = mlo(irow);  mi_up = 2.0 * mi_lo;

                    bndry.mi_lo = mi_lo;   bndry.mi_up = mi_up;
                    bndry.mj_lo = mj_lo;   bndry.mj_up = mj_up;
                    bndry.mjj   = [];      bndry.rjj   = [];
                    bndry.rvjj  = [];

                    b3(irow, jcol) = quadl(@(x) obj.integr3a(x, int_param, bndry), ...
                                           mi_lo, mi_up) / (mi_lo * mj_lo);
                end
            end

            % -------------------- Case 2: b2 (gain: i < j) --------------
            b2 = beta_init;
            warning('off'); % match legacy quiet integration
            for jcol = 2:n_sections
                for irow = 1:(jcol - 1)
                    mj_lo = mlo(jcol);  mj_up = 2.0 * mj_lo;
                    mi_lo = mlo(irow);  mi_up = 2.0 * mi_lo;

                    bndry.mi_lo = mi_lo;   bndry.mi_up = mi_up;
                    bndry.mj_lo = mj_lo;   bndry.mj_up = mj_up;
                    bndry.mjj   = [];      bndry.rjj   = [];
                    bndry.rvjj  = [];

                    b2(irow, jcol) = quadl(@(x) obj.integr2a(x, int_param, bndry), ...
                                           mi_lo, mi_up) / (mi_lo * mj_lo);
                end
            end
            warning('on');

            % -------------------- Case 1: b1 (gain: (j-1) & i < j) ------
            b1 = beta_init;
            for jcol = 2:n_sections
                for irow = 1:(jcol - 1)
                    mj_lo = mlo(jcol - 1);  mj_up = 2.0 * mj_lo;
                    mi_lo = mlo(irow);      mi_up = 2.0 * mi_lo;

                    bndry.mi_lo = mi_lo;   bndry.mi_up = mi_up;
                    bndry.mj_lo = mj_lo;   bndry.mj_up = mj_up;
                    bndry.mjj   = [];      bndry.rjj   = [];
                    bndry.rvjj  = [];

                    b1(irow, jcol) = quadl(@(x) obj.integr1a(x, int_param, bndry), ...
                                           mi_lo, mi_up) / (mi_lo * mj_lo);
                end
            end
            % super-diagonal double-count correction
            b1 = b1 - 0.5 * diag(diag(b1, 1), 1);

            % -------------------- Package into BetaMatrices -------------
            betas      = BetaMatrices();
            betas.b1   = b1;
            betas.b2   = b2;
            betas.b3   = b3;
            betas.b4   = b4;
            betas.b5   = b5;

            % >>> Consistent convention for RHS
            betas.b25  = b2 + b5;
        end

        % =================================================================
        %  INTEGRATION HELPERS (unchanged)
        % =================================================================

        function x = integr5a(obj, mj, param, bndry)
            nj  = length(mj);
            x   = 0 * mj;
            rj  = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj  = mj(iv);
                bndry.rjj  = rj(iv);
                bndry.rvjj = rvj(iv);
                x(iv) = quadl(@(y) obj.integr5b(y, param, bndry), ...
                              bndry.mj_lo, bndry.mj_up);
            end
            x = x ./ mj;
        end

        function yint = integr5b(~, mi, param, bndry)
            ri  = param.amfrac * mi.^(param.bmfrac);
            ni  = length(mi);
            rj  = bndry.rjj  * ones(1, ni);
            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
        end

        function x = integr4a(obj, mj, param, bndry)
            nj  = length(mj);
            x   = 0 * mj;
            rj  = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj  = mj(iv);
                bndry.rjj  = rj(iv);
                bndry.rvjj = rvj(iv);
                x(iv) = quadl(@(y) obj.integr4b(y, param, bndry), ...
                              bndry.mj_lo, bndry.mj_up);
            end
        end

        function yint = integr4b(~, mi, param, bndry)
            ri  = param.amfrac * mi.^(param.bmfrac);
            ni  = length(mi);
            rj  = bndry.rjj  * ones(1, ni);
            mj  = bndry.mjj  * ones(1, ni);
            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
            yint = (mi + mj) ./ mi ./ mj .* yint;
        end

        function x = integr3a(obj, mj, param, bndry)
            nj  = length(mj);
            x   = 0 * mj;
            rj  = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj  = mj(iv);
                bndry.rjj  = rj(iv);
                bndry.rvjj = rvj(iv);
                x(iv) = quadl(@(y) obj.integr3b(y, param, bndry), ...
                              bndry.mj_up - bndry.mjj, bndry.mj_up);
            end
            x = x ./ mj;
        end

        function yint = integr3b(~, mi, param, bndry)
            ri  = param.amfrac * mi.^(param.bmfrac);
            ni  = length(mi);
            rj  = bndry.rjj  * ones(1, ni);
            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
        end

        function x = integr2a(obj, mj, param, bndry)
            nj  = length(mj);
            x   = 0 * mj;
            rj  = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj  = mj(iv);
                bndry.rjj  = rj(iv);
                bndry.rvjj = rvj(iv);
                x(iv) = quadl(@(y) obj.integr2b(y, param, bndry), ...
                              bndry.mj_lo, bndry.mj_up - bndry.mjj);
            end
        end

        function yint = integr2b(~, mi, param, bndry)
            ri  = param.amfrac * mi.^(param.bmfrac);
            ni  = length(mi);
            rj  = bndry.rjj  * ones(1, ni);
            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
            yint = yint ./ mi;
        end

        function x = integr1a(obj, mj, param, bndry)
            nj  = length(mj);
            x   = 0 * mj;
            rj  = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj  = mj(iv);
                bndry.rjj  = rj(iv);
                bndry.rvjj = rvj(iv);
                mlow = max([bndry.mj_up - bndry.mjj, bndry.mj_lo]);
                x(iv) = quadl(@(y) obj.integr1b(y, param, bndry), ...
                              mlow, bndry.mj_up);
            end
        end

        function yint = integr1b(~, mi, param, bndry)
            ri  = param.amfrac * mi.^(param.bmfrac);
            ni  = length(mi);
            rj  = bndry.rjj  * ones(1, ni);
            mj  = bndry.mjj  * ones(1, ni);
            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
            yint = yint .* (mi + mj) ./ mi ./ mj;
        end
    end
end