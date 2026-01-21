classdef DiagnosticsSuite
%DIAGNOSTICSSUITE  Standard diagnostics for column coagulation runs.
%
% Part 1:
%   (1) Column inventory + budget curves
%   (2) Bottom export split into size classes (500/2000 um)
%   (3) Size spectra at selected depths
%
% Part 2:
%   (4) "Detective plots": dQ/dt vs size by-process at selected depths
%
% Notes:
% - Prefers sim.rhs.decomposeTerms(t,v) when available for consistent splitting.
% - Column integration assumes dv is a tendency in concentration units
%   consistent with inventory-density (biovolume per water volume).
% - Epsilon time series supports cfg.eps_fun OR cfg.epsilon_time+epsilon_series
%   OR cfg.epsilon_const.

    methods (Static)

        % ==========================================================
        % Public runners
        % ==========================================================
        function runPart1(simOrOut, outdir)
            if nargin < 2 || isempty(outdir)
                outdir = fullfile(pwd, 'diagnostics_part1');
            end
            if ~exist(outdir, 'dir'); mkdir(outdir); end

            [sim, out, cfg, grd] = DiagnosticsSuite.unpackSimOrOut(simOrOut);
            out = DiagnosticsSuite.ensureOutputData(out, grd, cfg);

            DiagnosticsSuite.plotMassBalanceCurves(sim, out, cfg, grd, outdir);
            DiagnosticsSuite.plotExportSizeClasses(out, cfg, grd, outdir);
            DiagnosticsSuite.plotSizeSpectra(out, cfg, grd, outdir);

            fprintf('\n[Part 1 DONE] Saved diagnostics to: %s\n', outdir);
        end


        function runPart2(sim, outdir)
            if nargin < 2 || isempty(outdir)
                outdir = fullfile(pwd, 'diagnostics_part2');
            end
            if ~exist(outdir, 'dir'); mkdir(outdir); end

            if ~isa(sim,'CoagulationSimulation')
                error('DiagnosticsSuite.runPart2 expects a CoagulationSimulation object.');
            end

            out = sim.result;
            cfg = sim.config;
            grd = sim.grid;

            out = DiagnosticsSuite.ensureOutputData(out, grd, cfg);

            DiagnosticsSuite.plotTendenciesByProcess(sim, out, cfg, grd, outdir);
            fprintf('\n[Part 2 DONE] Saved diagnostics to: %s\n', outdir);
        end


        % ==========================================================
        % (1) Mass balance curves (column-integrated)
        % ==========================================================
        function plotMassBalanceCurves(sim, out, cfg, grd, outdir)
            DiagnosticsSuite.assertColumnMode(cfg);

            t  = out.time(:);
            od = out.output_data;

            Ns    = cfg.n_sections;
            Nz    = cfg.getNumLayers();
            dz_cm = cfg.dz * 100;

            % Inventory [cm^3 m^-2]
            M = DiagnosticsSuite.pickVectorField(od, ...
                {'column_total_inventory_cm3m2','total_mass'}, ...
                'inventory');

            % Bottom export [cm^3 m^-2 d^-1] (+ down)
            F_export = DiagnosticsSuite.pickVectorField(od, ...
                {'bottom_total_flux_cm3m2d','total_flux'}, ...
                'bottom export flux');

            % FD reference
            dMdt_fd = gradient(M, t);

            % Growth operator (if present)
            G = [];
            if isfield(out,'operators') && isfield(out.operators,'growth') && ~isempty(out.operators.growth)
                G = out.operators.growth;
            end

            have_sim = ~isempty(sim) && isa(sim,'CoagulationSimulation') && ismethod(sim.rhs,'decomposeTerms');

            growth_rate = zeros(size(t));
            pp_rate     = zeros(size(t));
            coag_rate   = zeros(size(t));
            disagg_rate = zeros(size(t));
            nonlin_rate = zeros(size(t));

            for it = 1:numel(t)
                v = out.concentrations(it,:).';

                % Growth (if operator exists)
                if ~isempty(G)
                    dv_g = G * v;

                    % allow 0-D growth operator (Ns) applied only at surface layer
                    if numel(dv_g) == Ns
                        dv_full = zeros(Ns*Nz,1);
                        dv_full(1:Ns) = dv_g;
                        dv_g = dv_full;
                    end

                    growth_rate(it) = DiagnosticsSuite.integrateColumnRate(dv_g, Ns, Nz, dz_cm);
                end

                % PP/coag/disagg via decomposition (preferred)
                if have_sim
                    terms = sim.rhs.decomposeTerms(t(it), v);

                    pp_rate(it)     = DiagnosticsSuite.integrateColumnRate(terms.dv_pp,     Ns, Nz, dz_cm);
                    coag_rate(it)   = DiagnosticsSuite.integrateColumnRate(terms.dv_coag,   Ns, Nz, dz_cm);
                    disagg_rate(it) = DiagnosticsSuite.integrateColumnRate(terms.dv_disagg, Ns, Nz, dz_cm);
                    nonlin_rate(it) = coag_rate(it) + disagg_rate(it);
                end
            end

            % Residual sanity:
            % dM/dt ≈ growth + PP - export + nonlinear
            residual = dMdt_fd - growth_rate - pp_rate + F_export - nonlin_rate;

            % cumulative check
            cumInput  = cumtrapz(t, (growth_rate + pp_rate));
            cumExport = cumtrapz(t, F_export);
            cumNonlin = cumtrapz(t, nonlin_rate);
            dM        = M - M(1);

            fig = figure('Color','w','Position',[80 80 1080 900], 'Name','Part1 Mass Balance');
            DiagnosticsSuite.killToolbar(fig);

            subplot(3,1,1);
            plot(t, M, 'LineWidth',2);
            xlabel('Time [d]'); ylabel('Inventory [cm^3 m^{-2}]');
            title('Column Inventory');
            grid on; box on;

            subplot(3,1,2);
            plot(t, dMdt_fd, 'LineWidth',2); hold on;
            plot(t, growth_rate, '--', 'LineWidth',2);
            plot(t, pp_rate,     '--', 'LineWidth',2);
            plot(t, F_export,    '--', 'LineWidth',2);
            plot(t, nonlin_rate, ':',  'LineWidth',2);
            plot(t, residual,    '-.', 'LineWidth',2);

            legend('dM/dt (FD)','Growth/Input','PP source','Export (bottom)', ...
                   'Nonlinear (coag+disagg)','Residual', 'Location','best');

            xlabel('Time [d]'); ylabel('[cm^3 m^{-2} d^{-1}]');
            title('Rates: dM/dt ≈ Growth + PP − Export + Nonlinear');
            grid on; box on;

            subplot(3,1,3);
            plot(t, dM, 'LineWidth',2); hold on;
            plot(t, cumInput,  '--', 'LineWidth',2);
            plot(t, cumExport, '--', 'LineWidth',2);
            plot(t, cumNonlin, ':',  'LineWidth',2);
            plot(t, (cumInput - cumExport + cumNonlin), '-.', 'LineWidth',2);

            xlabel('Time [d]'); ylabel('[cm^3 m^{-2}]');
            legend('\DeltaM','\int (Growth+PP)','\int Export','\int Nonlinear', ...
                   '\int(Growth-Export+Nonlin)', 'Location','best');
            title('Cumulative Check');
            grid on; box on;

            DiagnosticsSuite.saveFig(fig, fullfile(outdir, 'part1_mass_balance_curves.png'));
            close(fig);
        end


        % ==========================================================
        % (2) Bottom export split by size classes
        % ==========================================================
        function plotExportSizeClasses(out, cfg, grd, outdir)
            t  = out.time(:);
            od = out.output_data;

            % Expect Fbin: Nt x Ns (bottom export per size bin)
            Fbin = DiagnosticsSuite.pickMatrixField(od, {'bottom_fluxsect_cm3m2d'}, 'bottom_fluxsect_cm3m2d');

            % Fix orientation if user stored as Ns x Nt
            if size(Fbin,1) ~= numel(t) && size(Fbin,2) == numel(t)
                Fbin = Fbin.'; % Nt x Ns
            end

            % Now enforce Nt x Ns
            Nt = numel(t);
            if size(Fbin,1) ~= Nt
                error('DiagnosticsSuite: bottom_fluxsect_cm3m2d must have Nt rows (Nt=%d). Got %dx%d.', ...
                    Nt, size(Fbin,1), size(Fbin,2));
            end

            Ns = cfg.n_sections;
            if size(Fbin,2) ~= Ns
                error(['DiagnosticsSuite: bottom_fluxsect_cm3m2d must have Ns columns (Ns=%d). ' ...
                       'Got %dx%d. Check OutputGenerator output shape.'], ...
                       Ns, size(Fbin,1), size(Fbin,2));
            end

            Ftot = sum(Fbin, 2);

            [D_um, D_label] = DiagnosticsSuite.getDiametersUm(od, grd);
            if numel(D_um) ~= Ns
                error('DiagnosticsSuite: diameter vector length (%d) must equal Ns (%d).', numel(D_um), Ns);
            end

            % Adrian bins: <500, 500–2000, >2000
            idx_small  = (D_um < 500);
            idx_medium = (D_um >= 500) & (D_um < 2000);
            idx_large  = (D_um >= 2000);

            F_small  = sum(Fbin(:, idx_small),  2);
            F_medium = sum(Fbin(:, idx_medium), 2);
            F_large  = sum(Fbin(:, idx_large),  2);

            Ftot_safe   = max(Ftot, 1e-30);
            frac_small  = F_small  ./ Ftot_safe;
            frac_medium = F_medium ./ Ftot_safe;
            frac_large  = F_large  ./ Ftot_safe;

            % steady scaling based on last 20%
            idx_ref = max(1, numel(t)-round(0.2*numel(t))):numel(t);
            F_ss = median(Ftot(idx_ref));
            if ~isfinite(F_ss) || abs(F_ss) < 1e-30; F_ss = 1e-30; end

            fig = figure('Color','w','Position',[120 80 1080 920], 'Name','Part1 Export Size Classes');
            DiagnosticsSuite.killToolbar(fig);

            subplot(3,1,1);
            eps_t = DiagnosticsSuite.getEpsilonTimeseries(cfg, t);
            if ~isempty(eps_t)
                semilogy(t, eps_t, 'k-', 'LineWidth',1.8);
                ylabel('\epsilon'); title('Turbulence Forcing');
            else
                plot(t, 0*t, 'k-'); ylim([-1 1]);
                ylabel('n/a'); title('Turbulence Forcing (disabled / not provided)');
            end
            grid on; box on;

            subplot(3,1,2);
            plot(t, Ftot/F_ss, 'k-', 'LineWidth',2); hold on;
            plot(t, F_small/F_ss,  '--', 'LineWidth',2);
            plot(t, F_medium/F_ss, '--', 'LineWidth',2);
            plot(t, F_large/F_ss,  '--', 'LineWidth',2);
            ylabel('Flux / steady');
            legend('Total','Small <500','Med 500–2000','Large >2000', 'Location','best');
            title(sprintf('Relative Bottom Export Flux (%s) | maxD=%.0f \\mum', D_label, max(D_um)));
            grid on; box on;

            subplot(3,1,3);
            plot(t, frac_large,  'LineWidth',2); hold on;
            plot(t, frac_medium, 'LineWidth',2);
            plot(t, frac_small,  'LineWidth',2);
            ylim([0 1]);
            ylabel('Fraction'); xlabel('Time [d]');
            legend('Large >2000','Med 500–2000','Small <500', 'Location','best');
            title('Bottom Export Fraction by Size Class');
            grid on; box on;

            DiagnosticsSuite.saveFig(fig, fullfile(outdir, 'part1_export_size_classes_500_2000.png'));
            close(fig);
        end


        % ==========================================================
        % (3) Size spectra: log N vs log D
        % ==========================================================
        function plotSizeSpectra(out, cfg, grd, outdir)
            DiagnosticsSuite.assertColumnMode(cfg);

            t  = out.time(:);
            od = out.output_data;

            Ns = cfg.n_sections;

            Y3 = DiagnosticsSuite.pickArrayField(od, {'Y3'}, 'Y3 (depth-resolved spectra)');
            [D_um, D_label] = DiagnosticsSuite.getDiametersUm(od, grd);

            v_lower = grd.v_lower(:);
            dwidth  = grd.dwidth(:);

            zc = DiagnosticsSuite.getDepthCenters(cfg, grd);

            depths_m = DiagnosticsSuite.pickDepths(zc, cfg);
            times_d  = DiagnosticsSuite.pickTimes(t, [0 5 10 18 19 20 21 22 23 25 30]);

            for id = 1:numel(depths_m)
                [~, kz] = min(abs(zc - depths_m(id)));

                fig = figure('Color','w','Position',[120 120 1040 720], ...
                    'Name',sprintf('Part1 SizeSpectra z=%.1fm', zc(kz)));
                DiagnosticsSuite.killToolbar(fig);

                hold on;
                leg = cell(numel(times_d),1);

                for jt = 1:numel(times_d)
                    tt = times_d(jt);
                    [~, it] = min(abs(t - tt));

                    vvol = squeeze(Y3(it, kz, :));  % Ns x 1
                    if numel(vvol) ~= Ns
                        error('DiagnosticsSuite: Y3 slice is not Ns-long at (it=%d,kz=%d).', it, kz);
                    end

                    nspec_v = vvol(:) ./ max(1.5*v_lower, eps) ./ max(dwidth, eps);
                    nspec_v = max(nspec_v, 1e-30);

                    plot(D_um, nspec_v, 'LineWidth',1.6);
                    leg{jt} = sprintf('t=%.2f d (snap %.2f)', tt, t(it));
                end

                set(gca,'XScale','log','YScale','log');
                xlabel(sprintf('%s [\\mum]', D_label));
                ylabel('Number spectrum (legacy style)');
                title(sprintf('Size spectrum at z \\approx %.1f m', zc(kz)));
                legend(leg, 'Location','eastoutside');
                grid on; box on;

                DiagnosticsSuite.saveFig(fig, fullfile(outdir, sprintf('part1_sizespectra_z_%dm.png', round(zc(kz)))));
                close(fig);
            end
        end


        % ==========================================================
        % (4) Part 2: dQ/dt vs size by-process
        % ==========================================================
        function plotTendenciesByProcess(sim, out, cfg, grd, outdir)
            DiagnosticsSuite.assertColumnMode(cfg);

            if ~ismethod(sim.rhs,'decomposeTerms')
                error('plotTendenciesByProcess requires sim.rhs.decomposeTerms().');
            end

            t  = out.time(:);
            od = out.output_data;

            Ns = cfg.n_sections;
            Nz = cfg.getNumLayers();

            [D_um, D_label] = DiagnosticsSuite.getDiametersUm(od, grd);
            zc = DiagnosticsSuite.getDepthCenters(cfg, grd);

            depths_m = DiagnosticsSuite.pickDepths(zc, cfg);
            times_d  = DiagnosticsSuite.pickTimes(t, [0 5 10 18 19 20 21 22 23 25]);

            titles = {'Total dQ/dt','Linear (growth-sinking)','Coagulation','PP source','Disaggregation'};
            fields = {'dv_tot','dv_lin','dv_coag','dv_pp','dv_disagg'};

            for id = 1:numel(depths_m)
                [~, kz] = min(abs(zc - depths_m(id)));

                fig = figure('Color','w','Position',[120 120 1260 820], ...
                    'Name',sprintf('Part2 Tendencies z=%.1fm', zc(kz)));
                DiagnosticsSuite.killToolbar(fig);

                for p = 1:5
                    subplot(2,3,p); hold on;

                    for jt = 1:numel(times_d)
                        tt = times_d(jt);
                        [~, it] = min(abs(t - tt));

                        v = out.concentrations(it,:).';
                        terms = sim.rhs.decomposeTerms(t(it), v);

                        dv = terms.(fields{p});
                        dv_layer = reshape(dv(:), [Ns, Nz]);
                        y = dv_layer(:, kz);

                        plot(D_um, y, 'LineWidth',1.4, 'DisplayName', sprintf('t=%.2f', tt));
                    end

                    set(gca,'XScale','log');
                    yline(0,'k-');
                    xlabel(sprintf('%s [\\mum]', D_label));
                    ylabel('dQ/dt [per day]');
                    title(titles{p});
                    grid on; box on;

                    if p == 1
                        legend('Location','best');
                    end
                end

                subplot(2,3,6);
                axis off;
                text(0.05, 0.70, sprintf('z \\approx %.1f m', zc(kz)), 'FontSize', 12);
                text(0.05, 0.50, sprintf('times: %s', strjoin(arrayfun(@(x) sprintf('%.1f',x), times_d, 'uni',0), ', ')));

                DiagnosticsSuite.saveFig(fig, fullfile(outdir, sprintf('part2_tendencies_z_%dm.png', round(zc(kz)))));
                close(fig);
            end
        end


        % ==========================================================
        % Helpers
        % ==========================================================
        function [sim, out, cfg, grd] = unpackSimOrOut(simOrOut)
            sim = [];
            if isa(simOrOut, 'CoagulationSimulation')
                sim = simOrOut;
                out = sim.result;
                cfg = sim.config;
                grd = sim.grid;
                return;
            end

            out = simOrOut;

            cfg = [];
            if isfield(out,'operators') && isfield(out.operators,'config')
                cfg = out.operators.config;
            end
            if isempty(cfg)
                error('DiagnosticsSuite: cannot find config. Pass the sim object instead.');
            end
            grd = DerivedGrid(cfg);
        end


        function out = ensureOutputData(out, grd, cfg)
            if ~isfield(out,'output_data') || isempty(out.output_data)
                out.output_data = OutputGenerator.spectraAndFluxes(out.time, out.concentrations, grd, cfg);
            end
        end


        function assertColumnMode(cfg)
            if ~isprop(cfg,'use_column') || ~cfg.use_column
                error('DiagnosticsSuite expects column mode (cfg.use_column=true).');
            end
        end


        % ---- Field pickers (IMPORTANT: do NOT flatten matrices!) ----
        function v = pickVectorField(S, names, what)
            v = [];
            for k = 1:numel(names)
                if isfield(S, names{k}) && ~isempty(S.(names{k}))
                    v = S.(names{k});
                    break;
                end
            end
            if isempty(v)
                error('DiagnosticsSuite: missing %s (%s).', what, strjoin(names, ', '));
            end
            if ~isvector(v)
                error('DiagnosticsSuite: %s must be a vector. Got %dx%d.', what, size(v,1), size(v,2));
            end
            v = v(:);
        end


        function A = pickMatrixField(S, names, what)
            A = [];
            for k = 1:numel(names)
                if isfield(S, names{k}) && ~isempty(S.(names{k}))
                    A = S.(names{k});
                    break;
                end
            end
            if isempty(A)
                error('DiagnosticsSuite: missing %s (%s).', what, strjoin(names, ', '));
            end
            if ~ismatrix(A)
                error('DiagnosticsSuite: %s must be a 2-D matrix. Got ndims=%d.', what, ndims(A));
            end
        end


        function X = pickArrayField(S, names, what)
            X = [];
            for k = 1:numel(names)
                if isfield(S, names{k}) && ~isempty(S.(names{k}))
                    X = S.(names{k});
                    break;
                end
            end
            if isempty(X)
                error('DiagnosticsSuite: missing %s (%s).', what, strjoin(names, ', '));
            end
        end


        function [D_um, D_label] = getDiametersUm(od, grd)
            if isfield(od,'diam_i') && ~isempty(od.diam_i)
                D_cm = od.diam_i(:);
                D_label = 'D_i';
            else
                D_cm = grd.getVolumeDiameters();
                D_label = 'D_v';
            end
            D_um = D_cm * 1e4;
        end


        function zc = getDepthCenters(cfg, grd)
            if isprop(grd,'z_centers') && ~isempty(grd.z_centers)
                zc = grd.z_centers(:);
            else
                zc = cfg.getZ(:);
            end
        end


        function depths_m = pickDepths(zc, cfg)
            zmax = max(zc);
            if zmax >= 100
                depths_m = [10, 50, 100];
            else
                z1 = min(10, zmax);
                z2 = min(0.75*zmax, zmax);
                z3 = max(zmax - cfg.dz/2, z1);
                depths_m = unique([z1, z2, z3], 'stable');
            end
        end


        function times_d = pickTimes(t, candidates)
            tmin = min(t); tmax = max(t);
            candidates = candidates(candidates >= tmin & candidates <= tmax);
            if isempty(candidates)
                candidates = unique([tmin, median(t), tmax]);
            end
            times_d = unique(candidates(:).');
        end


        function rate_cm3m2d = integrateColumnRate(dv, Ns, Nz, dz_cm)
            dv2 = reshape(dv(:), [Ns, Nz]);     % Ns x Nz
            layer_sum = sum(dv2, 1);            % 1 x Nz
            rate_cm3cm2d = sum(layer_sum) * dz_cm;
            rate_cm3m2d  = rate_cm3cm2d * 1e4;  % cm^2 -> m^2
        end


        function eps_t = getEpsilonTimeseries(cfg, t)
            eps_t = [];

            % (1) eps_fun(t,z)
            if isprop(cfg,'eps_fun') && ~isempty(cfg.eps_fun) && isa(cfg.eps_fun,'function_handle')
                z0 = 0;
                eps_t = arrayfun(@(tt) cfg.eps_fun(tt, z0), t(:));
                eps_t = max(eps_t, 1e-30);
                return;
            end

            % (2) epsilon_time + epsilon_series
            if isprop(cfg,'epsilon_time') && isprop(cfg,'epsilon_series') && ...
                    ~isempty(cfg.epsilon_time) && ~isempty(cfg.epsilon_series)

                tt = cfg.epsilon_time(:);
                E  = cfg.epsilon_series;

                if isvector(E)
                    ee = E(:);
                else
                    ee = E(:,1); % assume surface
                end

                eps_t = interp1(tt, ee, t(:), 'linear', 'extrap');
                eps_t = max(eps_t, 1e-30);
                return;
            end

            % (3) epsilon_const
            if isprop(cfg,'epsilon_const') && ~isempty(cfg.epsilon_const)
                eps_t = cfg.epsilon_const * ones(size(t(:)));
                eps_t = max(eps_t, 1e-30);
                return;
            end
        end


        function killToolbar(fig)
            try, fig.ToolBar = 'none'; catch, end
            try, set(fig, 'MenuBar','none'); catch, end
        end


        function saveFig(fig, filename)
            try
                exportgraphics(fig, filename, 'Resolution', 200);
            catch
                saveas(fig, filename);
            end
        end

    end
end