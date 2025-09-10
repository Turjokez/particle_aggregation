classdef BetaAssembler < handle
    %BETAASSEMBLER Computes sectionally integrated coagulation kernel matrices
    properties
        config;         % SimulationConfig object
        grid;           % DerivedGrid object
        kernels;        % KernelLibrary reference
    end

    methods
        function obj = BetaAssembler(config, grid)
            obj.config = config;
            obj.grid = grid;
            obj.kernels = KernelLibrary();
        end

        function betas = computeFor(obj, kernelName)
            original_kernel = obj.config.kernel;
            obj.config.kernel = kernelName;
            betas = obj.computeBetaMatrices();
            obj.config.kernel = original_kernel;
        end

        function betas = combineAndScale(obj, b_brown, b_shear, b_ds)
            % Combine and scale different kernel contributions
            % Scale Brownian by conBr*day_to_sec, shear by gamma*day_to_sec, DS by setcon*day_to_sec

            if ~isempty(b_brown) && ~isempty(b_brown.b1)
                b_brown_scaled = obj.scaleBetas(b_brown, obj.grid.conBr * obj.config.day_to_sec);
            else
                b_brown_scaled = BetaMatrices();
            end

            if ~isempty(b_shear) && ~isempty(b_shear.b1)
                b_shear_scaled = obj.scaleBetas(b_shear, obj.config.gamma * obj.config.day_to_sec);
            else
                b_shear_scaled = BetaMatrices();
            end

            if ~isempty(b_ds) && ~isempty(b_ds.b1)
                b_ds_scaled = obj.scaleBetas(b_ds, obj.grid.setcon * obj.config.day_to_sec);
            else
                b_ds_scaled = BetaMatrices();
            end

            % Combine all contributions
            betas = BetaMatrices();
            if ~isempty(b_brown_scaled.b1) || ~isempty(b_shear_scaled.b1) || ~isempty(b_ds_scaled.b1)
                % Initialize with zeros if any matrix is empty
                if isempty(b_brown_scaled.b1), b_brown_scaled = BetaMatrices(); end
                if isempty(b_shear_scaled.b1), b_shear_scaled = BetaMatrices(); end
                if isempty(b_ds_scaled.b1), b_ds_scaled = BetaMatrices(); end

                betas.b1 = b_brown_scaled.b1 + b_shear_scaled.b1 + b_ds_scaled.b1;
                betas.b2 = b_brown_scaled.b2 + b_shear_scaled.b2 + b_ds_scaled.b2;
                betas.b3 = b_brown_scaled.b3 + b_shear_scaled.b3 + b_ds_scaled.b3;
                betas.b4 = b_brown_scaled.b4 + b_shear_scaled.b4 + b_ds_scaled.b4;
                betas.b5 = b_brown_scaled.b5 + b_shear_scaled.b5 + b_ds_scaled.b5;
                betas.b25 = betas.b2 - betas.b3 - betas.b4 - betas.b5;
            else
                % If no matrices provided, create empty ones
                betas = BetaMatrices();
            end
        end

        function betas_scaled = scaleBetas(obj, betas, scale_factor)
            % Scale all beta matrices by a factor
            betas_scaled = BetaMatrices();
            betas_scaled.b1 = betas.b1 * scale_factor;
            betas_scaled.b2 = betas.b2 * scale_factor;
            betas_scaled.b3 = betas.b3 * scale_factor;
            betas_scaled.b4 = betas.b4 * scale_factor;
            betas_scaled.b5 = betas.b5 * scale_factor;
            betas_scaled.b25 = betas.b25 * scale_factor;
        end
    end

    methods (Access = private)
        function betas = computeBetaMatrices(obj)
            % Complete implementation ported from CalcBetas.m
            n_sections = obj.config.n_sections;
            mlo = obj.grid.v_lower;

            % Initialize beta matrices
            beta_init = zeros(n_sections, n_sections);

            % Set up integration parameters with proper structure for kernel functions
            int_param.amfrac = obj.grid.amfrac;
            int_param.bmfrac = obj.grid.bmfrac;
            int_param.kernel = obj.config.kernel;
            int_param.r_to_rg = obj.config.r_to_rg;
            int_param.setcon = obj.grid.setcon;
            int_param.constants = obj.config;

            % Case 5: loss from jcol by collisions of jcol & irow > jcol
            b5 = beta_init;
            for jcol = 1:(n_sections - 1)
                for irow = (jcol + 1):n_sections
                    mj_lo = mlo(jcol);
                    mj_up = 2.0 * mj_lo;
                    mi_lo = mlo(irow);
                    mi_up = 2.0 * mi_lo;

                    bndry.mi_lo = mi_lo;
                    bndry.mi_up = mi_up;
                    bndry.mj_lo = mj_lo;
                    bndry.mj_up = mj_up;
                    bndry.mjj = [];
                    bndry.rjj = [];
                    bndry.rvj = [];

                    b5(irow, jcol) = quadl(@(x) obj.integr5a(x, int_param, bndry), mi_lo, mi_up) / (mi_lo * mj_lo);
                end
            end

            % Case 4: loss from jcol by collisions with itself
            b4 = beta_init;
            for jcol = 1:n_sections
                mj_lo = mlo(jcol);
                mj_up = 2.0 * mj_lo;
                mi_lo = mlo(jcol);
                mi_up = 2.0 * mi_lo;

                bndry.mi_lo = mi_lo;
                bndry.mi_up = mi_up;
                bndry.mj_lo = mj_lo;
                bndry.mj_up = mj_up;
                bndry.mjj = [];
                bndry.rjj = [];
                bndry.rvj = [];

                b4(jcol, jcol) = quadl(@(x) obj.integr4a(x, int_param, bndry), mi_lo, mi_up) / (mi_lo * mj_lo);
            end
            % Take account of double counting by dividing by 2
            b4 = b4 / 2;

            % Case 3: loss from jcol by collisions of jcol & irow < jcol
            b3 = beta_init;
            for jcol = 2:n_sections
                for irow = 1:(jcol - 1)
                    mj_lo = mlo(jcol);
                    mj_up = 2.0 * mj_lo;
                    mi_lo = mlo(irow);
                    mi_up = 2.0 * mi_lo;

                    bndry.mi_lo = mi_lo;
                    bndry.mi_up = mi_up;
                    bndry.mj_lo = mj_lo;
                    bndry.mj_up = mj_up;
                    bndry.mjj = [];
                    bndry.rjj = [];
                    bndry.rvj = [];

                    b3(irow, jcol) = quadl(@(x) obj.integr3a(x, int_param, bndry), mi_lo, mi_up) / (mi_lo * mj_lo);
                end
            end

            % Case 2: gain in jcol by collisions of jcol & irow < jcol
            b2 = beta_init;
            warning('off'); % Suppress warnings during integration
            for jcol = 2:n_sections
                for irow = 1:(jcol - 1)
                    mj_lo = mlo(jcol);
                    mj_up = 2.0 * mj_lo;
                    mi_lo = mlo(irow);
                    mi_up = 2.0 * mi_lo;

                    bndry.mi_lo = mi_lo;
                    bndry.mi_up = mi_up;
                    bndry.mj_lo = mj_lo;
                    bndry.mj_up = mj_up;
                    bndry.mjj = [];
                    bndry.rjj = [];
                    bndry.rvj = [];

                    b2(irow, jcol) = quadl(@(x) obj.integr2a(x, int_param, bndry), mi_lo, mi_up) / (mi_lo * mj_lo);
                end
            end
            warning('on');

            % Case 1: gain in jcol by collisions of (jcol-1) & irow < jcol
            b1 = beta_init;
            for jcol = 2:n_sections
                for irow = 1:(jcol - 1)
                    mj_lo = mlo(jcol - 1);
                    mj_up = 2.0 * mj_lo;
                    mi_lo = mlo(irow);
                    mi_up = 2.0 * mi_lo;

                    bndry.mi_lo = mi_lo;
                    bndry.mi_up = mi_up;
                    bndry.mj_lo = mj_lo;
                    bndry.mj_up = mj_up;
                    bndry.mjj = [];
                    bndry.rjj = [];
                    bndry.rvj = [];

                    b1(irow, jcol) = quadl(@(x) obj.integr1a(x, int_param, bndry), mi_lo, mi_up) / (mi_lo * mj_lo);
                end
            end

            % Take account of double counting on the super-diagonal
            b1 = b1 - 0.5 * diag(diag(b1, 1), 1);

            % Create BetaMatrices object
            betas = BetaMatrices();
            betas.b1 = b1;
            betas.b2 = b2;
            betas.b3 = b3;
            betas.b4 = b4;
            betas.b5 = b5;
            betas.b25 = b2 - b3 - b4 - b5;
        end

        % Integration functions for case 5
        function x = integr5a(obj, mj, param, bndry)
            nj = length(mj);
            x = 0 * mj;
            rj = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj = mj(iv);
                bndry.rjj = rj(iv);
                bndry.rvjj = rvj(iv);

                x(iv) = quadl(@(y) obj.integr5b(y, param, bndry), bndry.mj_lo, bndry.mj_up);
            end
            x = x ./ mj;
        end

        function yint = integr5b(obj, mi, param, bndry)
            ri = param.amfrac * mi.^(param.bmfrac);

            ni = length(mi);
            rj = bndry.rjj * ones(1, ni);
            mj = bndry.mjj * ones(1, ni);

            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            % Get the actual kernel function handle
            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
        end

        % Integration functions for case 4
        function x = integr4a(obj, mj, param, bndry)
            nj = length(mj);
            x = 0 * mj;
            rj = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj = mj(iv);
                bndry.rjj = rj(iv);
                bndry.rvjj = rvj(iv);

                x(iv) = quadl(@(y) obj.integr4b(y, param, bndry), bndry.mj_lo, bndry.mj_up);
            end
        end

        function yint = integr4b(obj, mi, param, bndry)
            ri = param.amfrac * mi.^(param.bmfrac);

            ni = length(mi);
            rj = bndry.rjj * ones(1, ni);
            mj = bndry.mjj * ones(1, ni);

            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            % Get the actual kernel function handle
            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
            yint = (mi + bndry.mjj) ./ mi ./ bndry.mjj .* yint;
        end

        % Integration functions for case 3
        function x = integr3a(obj, mj, param, bndry)
            nj = length(mj);
            x = 0 * mj;
            rj = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj = mj(iv);
                bndry.rjj = rj(iv);
                bndry.rvjj = rvj(iv);

                x(iv) = quadl(@(y) obj.integr3b(y, param, bndry), bndry.mj_up - bndry.mjj, bndry.mj_up);
            end
            x = x ./ mj;
        end

        function yint = integr3b(obj, mi, param, bndry)
            ri = param.amfrac * mi.^(param.bmfrac);

            ni = length(mi);
            rj = bndry.rjj * ones(1, ni);
            mj = bndry.mjj * ones(1, ni);

            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            % Get the actual kernel function handle
            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
        end

        % Integration functions for case 2
        function x = integr2a(obj, mj, param, bndry)
            nj = length(mj);
            x = 0 * mj;
            rj = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj = mj(iv);
                bndry.rjj = rj(iv);
                bndry.rvjj = rvj(iv);

                x(iv) = quadl(@(y) obj.integr2b(y, param, bndry), bndry.mj_lo, bndry.mj_up - bndry.mjj);
            end
        end

        function yint = integr2b(obj, mi, param, bndry)
            ri = param.amfrac * mi.^(param.bmfrac);

            ni = length(mi);
            rj = bndry.rjj * ones(1, ni);
            mj = bndry.mjj * ones(1, ni);

            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            % Get the actual kernel function handle
            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
            yint = yint ./ mi;
        end

        % Integration functions for case 1
        function x = integr1a(obj, mj, param, bndry)
            nj = length(mj);
            x = 0 * mj;
            rj = param.amfrac * mj.^param.bmfrac;
            rvj = (0.75/pi * mj).^(1.0/3.0);

            for iv = 1:nj
                bndry.mjj = mj(iv);
                bndry.rjj = rj(iv);
                bndry.rvjj = rvj(iv);
                mlow = max([bndry.mj_up - bndry.mjj, bndry.mj_lo]);

                x(iv) = quadl(@(y) obj.integr1b(y, param, bndry), mlow, bndry.mj_up);
            end
        end

        function yint = integr1b(obj, mi, param, bndry)
            ri = param.amfrac * mi.^(param.bmfrac);

            ni = length(mi);
            rj = bndry.rjj * ones(1, ni);
            mj = bndry.mjj * ones(1, ni);

            rvi = (0.75/pi * mi).^(1.0/3.0);
            rvj = bndry.rvjj * ones(1, ni);

            % Get the actual kernel function handle
            kernel_func = KernelLibrary.getKernel(param.kernel);
            yint = kernel_func([ri; rj], [rvi; rvj], param);
            yint = yint .* (mi + mj) ./ mi ./ mj;
        end
    end
end
