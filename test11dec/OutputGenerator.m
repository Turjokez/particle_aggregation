classdef OutputGenerator < handle
    %OUTPUTGENERATOR Generates outputs and visualizations for coagulation simulation

    methods (Static)
        function output_data = spectraAndFluxes(t, Y, grid, config)
            %SPECTRAANDFLUXES Compute particle spectra and fluxes

            n_times = length(t);
            n_sections = length(grid.v_lower);

            % ==========================================================
            % Column-aware branch
            % ==========================================================
            if isprop(config, 'use_column') && config.use_column
                Nz = config.getNumLayers();
                Ns = config.n_sections;

                % Safety check
                if size(Y, 2) ~= Ns * Nz
                    error('OutputGenerator: Column Y dimension mismatch. Expected %d columns, got %d.', ...
                        Ns*Nz, size(Y,2));
                end

                % ------------------------------------------------------
                % Reshape Y to Y3(it, z, sec) = N(z,sec)
                % ------------------------------------------------------
                Y3 = zeros(n_times, Nz, Ns);
                for it = 1:n_times
                    % State ordering assumed:
                    % v = [N(:,1); N(:,2); ...; N(:,Nz)]  where N is Ns x Nz
                    N2 = reshape(Y(it, :)', [Ns, Nz]);     % Ns x Nz (sec x z)
                    Y3(it, :, :) = permute(N2, [2 1]);    % Nz x Ns
                end

                % Radii/diameters (same as 0-D)
                r_i   = grid.getFractalRadii();
                r_v   = grid.getConservedRadii();
                diam_i = grid.getImageDiameters(config);
                diam_v = grid.getVolumeDiameters();

                % Settling velocities:
                % velocity(...) returns [cm/s] in your legacy usage
                set_vel_cms   = SettlingVelocityService.velocity(r_i, r_v, grid.setcon); % cm/s
                set_vel_cmday = set_vel_cms * config.day_to_sec;                          % cm/day

                % Keep legacy-style row vectors for output fields
                diam_i = diam_i(:)';   % 1 x Ns
                diam_v = diam_v(:)';   % 1 x Ns

                diam_i_mat = diam_i(ones(n_times, 1), :);
                diam_v_mat = diam_v(ones(n_times, 1), :);

                % ------------------------------------------------------
                % STRICT closure-safe scalings
                % ------------------------------------------------------
                dz_m  = config.dz;          % m
                dz_cm = dz_m * 100;         % cm

                % Inventory per area:
                %   M = sum_z N(z) * dz_cm  -> cm^3/cm^2
                % then multiply by 1e4 to convert to cm^3/m^2
                column_inventory_by_section_cm3m2 = zeros(n_times, Ns);
                column_total_inventory_cm3m2      = zeros(n_times, 1);

                % Bottom export flux per area:
                %   F = N_bottom * w_cmday -> cm^3/cm^2/day
                % then *1e4 -> cm^3/m^2/day
                bottom_fluxsect_cm3m2d   = zeros(n_times, Ns);
                bottom_total_flux_cm3m2d = zeros(n_times, 1);

                % snapshots
                N_bottom  = zeros(n_times, Ns);
                N_surface = zeros(n_times, Ns);

                % (optional) legacy fields kept for plotting compatibility
                nspec_v_bottom   = zeros(n_times, Ns);
                masspec_v_bottom = zeros(n_times, Ns);
                fluxspec_bottom  = zeros(n_times, Ns);  %#ok<NASGU> % legacy-ish

                for it = 1:n_times
                    Nlayer = squeeze(Y3(it, :, :));   % Nz x Ns

                    ysurf = Nlayer(1, :);             % 1 x Ns
                    ybot  = Nlayer(end, :);           % 1 x Ns

                    N_surface(it, :) = ysurf;
                    N_bottom(it, :)  = ybot;

                    % ---- closure-safe inventory ----
                    column_inventory_by_section_cm3m2(it, :) = (sum(Nlayer, 1) * dz_cm) * 1e4;
                    column_total_inventory_cm3m2(it) = sum(column_inventory_by_section_cm3m2(it, :));

                    % ---- closure-safe bottom flux ----
                    bottom_fluxsect_cm3m2d(it, :)  = (ybot .* set_vel_cmday(:)') * 1e4;
                    bottom_total_flux_cm3m2d(it)   = sum(bottom_fluxsect_cm3m2d(it, :));

                    % ---- legacy spectra on bottom layer (OK to keep) ----
                    nspec_v_bottom(it, :)   = ybot ./ (1.5 * grid.v_lower') ./ grid.dwidth';
                    masspec_v_bottom(it, :) = ybot ./ grid.dwidth';
                end

                % cumulative export
                export_cum_cm3m2 = cumtrapz(t(:), bottom_total_flux_cm3m2d(:));

                % image/volume conversion ratio (legacy)
                diaratio = (config.fr_dim/3) * diam_v_mat ./ diam_i_mat;
                nspec_i_bottom   = nspec_v_bottom   .* diaratio;
                masspec_i_bottom = masspec_v_bottom .* diaratio;

                % Assemble output
                output_data = struct();

                output_data.Y3 = Y3;
                output_data.Nz = Nz;
                output_data.Ns = Ns;
                output_data.z  = config.getZ();
                output_data.dz = config.dz;

                % Closure-safe core fields
                output_data.column_inventory_by_section_cm3m2 = column_inventory_by_section_cm3m2;
                output_data.column_total_inventory_cm3m2      = column_total_inventory_cm3m2;

                output_data.bottom_fluxsect_cm3m2d   = bottom_fluxsect_cm3m2d;
                output_data.bottom_total_flux_cm3m2d = bottom_total_flux_cm3m2d;
                output_data.export_cum_cm3m2         = export_cum_cm3m2;

                % For backward compatibility with your plotting + test1 code
                output_data.bottom_total_flux = bottom_total_flux_cm3m2d;
                output_data.total_flux        = bottom_total_flux_cm3m2d;
                output_data.total_mass        = column_total_inventory_cm3m2; % NOTE: now "inventory", not raw sum(Y)

                % snapshots
                output_data.N_bottom  = N_bottom;
                output_data.N_surface = N_surface;

                % spectra (bottom layer)
                output_data.nspec_v   = nspec_v_bottom;
                output_data.masspec_v = masspec_v_bottom;
                output_data.nspec_i   = nspec_i_bottom;
                output_data.masspec_i = masspec_i_bottom;

                % diameters + grid info
                output_data.diam_i = diam_i(:);
                output_data.diam_v = diam_v(:);
                output_data.diam_i_mat = diam_i_mat;
                output_data.diam_v_mat = diam_v_mat;

                output_data.set_vel_cmday = set_vel_cmday(:);     % cm/day
                output_data.v_lower = grid.v_lower;
                output_data.dwidth  = grid.dwidth;
                output_data.diaratio = diaratio;
                output_data.fr_dim = config.fr_dim;

                % Optional debug save
                try
                    diag = struct();
                    diag.t = t;
                    diag.Y = Y;
                    diag.Y3 = Y3;
                    diag.set_vel_cmday = set_vel_cmday;
                    diag.M_cm3m2 = column_total_inventory_cm3m2;
                    diag.F_cm3m2d = bottom_total_flux_cm3m2d;
                    diag.ExportCum_cm3m2 = export_cum_cm3m2;
                    save('plot_diag_oop.mat','diag');
                catch
                end

                return;
            end

            % ==========================================================
            % ORIGINAL 0-D slab code (unchanged below)
            % ==========================================================

            nspec_v = zeros(n_times, n_sections);
            masspec_v = nspec_v;
            fluxsect = nspec_v;
            fluxspec = nspec_v;

            r_i = grid.getFractalRadii();
            r_v = grid.getConservedRadii();
            diam_i = grid.getImageDiameters(config);
            diam_v = grid.getVolumeDiameters();

            set_vel = SettlingVelocityService.velocity(r_i, r_v, grid.setcon);
            set_vel = set_vel / 100 * config.day_to_sec;

            diam_i = diam_i';
            diam_v = diam_v';
            diam_i_mat = diam_i(ones(n_times, 1), :);
            diam_v_mat = diam_v(ones(n_times, 1), :);

            for jindx = 1:n_times
                yout = Y(jindx, :);

                nspec_v(jindx, :) = yout ./ (1.5 * grid.v_lower') ./ grid.dwidth';
                masspec_v(jindx, :) = yout ./ grid.dwidth';
                fluxsect(jindx, :) = yout .* set_vel' * 1e6;
                fluxspec(jindx, :) = masspec_v(jindx, :) .* set_vel' * 1e6;
            end

            total_flux = sum(fluxsect, 2);
            total_mass = sum(Y, 2);

            diaratio = (config.fr_dim/3) * diam_v_mat ./ diam_i_mat;
            nspec_i = nspec_v .* diaratio;
            masspec_i = masspec_v .* diaratio;
            fluxspec_i = fluxspec .* diaratio;

            output_data = struct();
            output_data.nspec_v = nspec_v;
            output_data.masspec_v = masspec_v;
            output_data.fluxsect = fluxsect;
            output_data.fluxspec = fluxspec;
            output_data.total_flux = total_flux;
            output_data.total_mass = total_mass;
            output_data.diam_i = diam_i;
            output_data.diam_v = diam_v;
            output_data.diam_i_mat = diam_i_mat;
            output_data.diam_v_mat = diam_v_mat;
            output_data.nspec_i = nspec_i;
            output_data.masspec_i = masspec_i;
            output_data.fluxspec_i = fluxspec_i;
            output_data.set_vel = set_vel;
            output_data.v_lower = grid.v_lower;
            output_data.dwidth = grid.dwidth;
            output_data.diaratio = diaratio;
            output_data.fr_dim = config.fr_dim;

            try
                diag = struct();
                diag.t = t;
                diag.Y = Y;
                diag.r_i = r_i;
                diag.r_v = r_v;
                diag.diam_i = diam_i;
                diag.diam_v = diam_v;
                diag.set_vel = set_vel;
                diag.v_lower = grid.v_lower;
                diag.dwidth = grid.dwidth;
                diag.nspec_v = nspec_v;
                diag.masspec_v = masspec_v;
                diag.fluxspec = fluxspec;
                diag.diaratio = (config.fr_dim/3) * diam_v_mat ./ diam_i_mat;
                diag.fluxspec_i = fluxspec_i;
                save('plot_diag_oop.mat','diag');
            catch
            end
        end

        % --- everything else unchanged (plotAll, plotSpectraAndFluxes, etc.) ---
        function plotAll(t, Y, output_data, gains, losses, config, sectional_gains, sectional_losses, betas)
            OutputGenerator.plotSpectraAndFluxes(t, Y, output_data);
            OutputGenerator.plotMassBalance(t, gains, losses);
            OutputGenerator.plotCoagVsSettRatio(t, losses);
            OutputGenerator.plot3DSurfaces(t, Y, sectional_gains, sectional_losses, output_data);
        end

        function plotCoagVsSettRatio(t, losses)
            figure(3);
            if isfield(losses, 'coag') && isfield(losses, 'sett')
                ratio = losses.coag ./ max(losses.sett, eps);
                plot(t, ratio);
                set(gca, 'FontName', 'Helvetica', 'FontSize', 14)
                xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 14)
                ylabel('(Coag Losses)/(Settling Losses)', 'FontName', 'Helvetica', 'FontSize', 14)
            end
        end

        function plotSpectraAndFluxes(t, Y, output_data)
            figure(1);

            subplot(2, 2, 1);
            ns_init = output_data.nspec_i(1, :);
            ns_final = output_data.nspec_i(end, :);
            loglog(output_data.diam_i, ns_init, 'b', ...
                output_data.diam_i, ns_final, 'r');
            xlabel('Particle diameter [cm]');
            ylabel('Number spectrum [# cm^{-4}]');
            axis tight;

            subplot(2, 2, 3);
            plot(output_data.diam_i_mat(:,2:end)', output_data.fluxspec_i(:,2:end)');
            xlabel('Particle image diameter [cm]');
            ylabel('Volume flux spectra [cm^2 m^{-2} d^{-1}]');
            axis tight;

            subplot(2, 4, 3);
            semilogy(t, Y, t, output_data.total_mass, '*--');
            xlabel('Time [d]');
            ylabel('Sectional concentration [vol/vol/sect]');

            subplot(2, 4, 7);
            if length(t) == size(output_data.fluxsect, 1) && length(t) == length(output_data.total_flux)
                plot(t, output_data.fluxsect, t, output_data.total_flux, '*--');
            else
                plot(t, output_data.fluxsect(1:length(t), :), t, output_data.total_flux(1:length(t)), '*--');
            end
            xlabel('Time [d]');
            ylabel('Sectional Flux [cm^3 m^{-2} d^{-1} sect^{-1}]');

            subplot(2, 4, 8);
            if length(t) == length(output_data.total_flux) && length(t) == length(output_data.total_mass)
                plot(t, output_data.total_flux ./ output_data.total_mass / 1e6);
            else
                plot(t, output_data.total_flux(1:length(t)) ./ output_data.total_mass(1:length(t)) / 1e6);
            end
            xlabel('Time [d]');
            ylabel('Average v [m d^{-1}]');
        end

        function plotMassBalance(t, gains, losses)
            figure(2);
            clf;
            if isfield(gains, 'growth') && isfield(losses, 'sett') && isfield(losses, 'coag')
                denominator = losses.sett + losses.coag;
                ratio = gains.growth ./ max(denominator, eps);
                plot(t, ratio);
                xlabel('Time [d]');
                ylabel('Gains/Losses');
                title('Total System Mass Balance');
            end
        end

        function plotConcentrationEvolution(t, Y, total_mass)
            figure(5);
            if size(Y, 1) == length(t) && size(Y, 2) > 0
                semilogy(t, Y);
                hold on;
                if length(total_mass) == length(t)
                    semilogy(t, total_mass, '*--', 'LineWidth', 2);
                end
                hold off;
                xlabel('Time [d]');
                ylabel('Sectional concentration [vol/vol/sect]');
                legend([arrayfun(@(i) sprintf('Section %d', i), 1:size(Y, 2), 'UniformOutput', false), {'Total Mass'}]);
            else
                warning('Dimension mismatch in plotConcentrationEvolution: t=%d, Y=%dx%d, total_mass=%d', ...
                    length(t), size(Y, 1), size(Y, 2), length(total_mass));
            end
        end

        function plot3DSurfaces(t, Y, gains, losses, output_data)
            n_sections = size(Y, 2);
            if length(t) ~= size(Y, 1) || n_sections == 0
                warning('Dimension mismatch in plot3DSurfaces: t=%d, Y=%dx%d', length(t), size(Y, 1), n_sections);
                return;
            end

            t_mat = t(:, ones(1, n_sections));
            if length(output_data.v_lower) == n_sections
                v_mat = output_data.v_lower';
                v_mat = v_mat(ones(length(t), 1), :);
            else
                warning('Dimension mismatch: v_lower length=%d, n_sections=%d', length(output_data.v_lower), n_sections);
                v_mat = ones(length(t), n_sections);
            end

            if size(gains, 1) == length(t) && size(gains, 2) == n_sections && ...
                    size(losses, 1) == length(t) && size(losses, 2) == n_sections

                ratio = gains ./ max(losses, eps);
                ratio(losses <= 0) = NaN;
                ratio(~isfinite(ratio)) = NaN;
                ratio(ratio < 0) = NaN;

                figure(4);
                surf(log10(v_mat), t_mat, ratio);
                xlabel('Volume');
                ylabel('Time [d]');
                zlabel('Gains/Losses');

                valid_ratios = ratio(~isnan(ratio) & ratio > 0);
                if ~isempty(valid_ratios)
                    zlim([0, max(valid_ratios)]);
                end
            end
        end

        function exportData(output_data, filename)
            if nargin < 2
                filename = sprintf('coagulation_output_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
            end
            save(filename, 'output_data');
            fprintf('Output data saved to: %s\n', filename);
        end
    end
end