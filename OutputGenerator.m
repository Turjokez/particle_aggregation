classdef OutputGenerator < handle
    %OUTPUTGENERATOR Generates outputs and visualizations for coagulation simulation

    methods (Static)
        function output_data = spectraAndFluxes(t, Y, grid, config)
            %SPECTRAANDFLUXES Compute particle spectra and fluxes
            % t = time vector
            % Y = concentration matrix
            % grid = DerivedGrid object
            % config = SimulationConfig object
            % Returns: struct with computed data

            n_times = length(t);
            n_sections = length(grid.v_lower);

            % Initialize output arrays
            nspec_v = zeros(n_times, n_sections);
            masspec_v = nspec_v;
            fluxsect = nspec_v;
            fluxspec = nspec_v;

            % Get radii and diameters
            r_i = grid.getFractalRadii();
            r_v = grid.getConservedRadii();
            diam_i = grid.getImageDiameters(config);
            diam_v = grid.getVolumeDiameters();

            % Calculate settling velocities
            set_vel = SettlingVelocityService.velocity(r_i, r_v, grid.setcon);
            set_vel = set_vel / 100 * config.day_to_sec;

            % Create time matrices for vectorized operations
            % diam_i and diam_v should be row vectors to match legacy
            diam_i = diam_i';
            diam_v = diam_v';
            diam_i_mat = diam_i(ones(n_times, 1), :);
            diam_v_mat = diam_v(ones(n_times, 1), :);

            % Compute spectra and fluxes for each time
            for jindx = 1:n_times
                yout = Y(jindx, :);

                nspec_v(jindx, :) = yout ./ (1.5 * grid.v_lower') ./ grid.dwidth';
                masspec_v(jindx, :) = yout ./ grid.dwidth';
                fluxsect(jindx, :) = yout .* set_vel' * 1e6;
                fluxspec(jindx, :) = masspec_v(jindx, :) .* set_vel' * 1e6;
            end

            % Total quantities
            total_flux = sum(fluxsect, 2);
            total_mass = sum(Y, 2);

            % Diameter ratio and image-based spectra
            diaratio = (config.fr_dim/3) * diam_v_mat ./ diam_i_mat;
            nspec_i = nspec_v .* diaratio;
            masspec_i = masspec_v .* diaratio;
            fluxspec_i = fluxspec .* diaratio;

            % Assemble output data
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

            % Diagnostics snapshot for cross-checking legacy vs OOP
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
                % Also store total gains/losses used for Fig 2 & 3 parity checks
                [tg, tl] = MassBalanceAnalyzer.total(Y, struct('sink_loss', diag(1).v_lower*0 + diag.set_vel(1))); %#ok<NASGU>
                save('plot_diag_oop.mat','diag');
            catch
            end
        end

        function plotAll(t, Y, output_data, gains, losses, config, sectional_gains, sectional_losses, betas)
            %PLOTALL Generate all standard plots
            % t = time vector
            % Y = concentration matrix
            % output_data = output from spectraAndFluxes
            % gains, losses = from MassBalanceAnalyzer
            % config = SimulationConfig object

            % Figure 1: Spectra and fluxes (match legacy axes/layout)
            OutputGenerator.plotSpectraAndFluxes(t, Y, output_data);

            % Figure 2: Mass balance (single ratio panel only, like legacy)
            OutputGenerator.plotMassBalance(t, gains, losses);

            % Figure 3: (Coag Losses)/(Settling Losses) vs time (legacy)
            OutputGenerator.plotCoagVsSettRatio(t, losses);

            % Figure 4: Legacy Figure 6 → sectional (Gains/Losses) surface
            % Render only the ratio surface in Figure 4
            OutputGenerator.plot3DSurfaces(t, Y, sectional_gains, sectional_losses, output_data);
        end

        function plotCoagVsSettRatio(t, losses)
            %PLOTCOAGVSSEtTRATIO Legacy Figure 3: (Coag Losses)/(Settling Losses)
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
            %PLOTSPECTRAANDFLUXES Plot particle spectra and fluxes
            % Match legacy layout exactly: mixed 2x2 and 2x4 subplot layout
            figure(1);

            % Number spectrum (subplot 2,2,1 - overlaid on 2x4 layout)
            subplot(2, 2, 1);
            % Use the precomputed nspec_i data that was calculated in spectraAndFluxes()
            % This matches the legacy calculation: nspec_i = nspec_v .* diaratio
            ns_init = output_data.nspec_i(1, :);
            ns_final = output_data.nspec_i(end, :);
            loglog(output_data.diam_i, ns_init, 'b', ...
                output_data.diam_i, ns_final, 'r');
            xlabel('Particle diameter [cm]');
            ylabel('Number spectrum [# cm^{-4}]');
            axis tight;

            % Flux spectra (subplot 2,2,3 - overlaid on 2x4 layout)
            subplot(2, 2, 3);
            % Use the same approach as legacy: plot(diam_i_mat(:,2:end)', fluxspec_i(:,2:end)')
            % This plots each time step as a separate line
            % diam_i_mat is n_times x n_sections, so (:,2:end)' is (n_sections-1) x n_times
            % fluxspec_i is n_times x n_sections, so (:,2:end)' is (n_sections-1) x n_times
            plot(output_data.diam_i_mat(:,2:end)', output_data.fluxspec_i(:,2:end)');
            xlabel('Particle image diameter [cm]');
            ylabel('Volume flux spectra [cm^2 m^{-2} d^{-1}]');
            axis tight;

            % Sectional concentration vs time (subplot 2,4,3 - this overwrites the 2,2,3 flux spectra)
            subplot(2, 4, 3);
            % Match legacy exactly: semilogy(t_out, spec, t_out, total_mass, '*--')
            semilogy(t, Y, t, output_data.total_mass, '*--');
            xlabel('Time [d]');
            ylabel('Sectional concentration [vol/vol/sect]');

            % Sectional flux over time with total (subplot 2,4,7)
            subplot(2, 4, 7);
            % Match legacy exactly: plot(t_out, fluxsect, t_out, total_flux, '*--')
            % Ensure dimensions match for plotting
            if length(t) == size(output_data.fluxsect, 1) && length(t) == length(output_data.total_flux)
                plot(t, output_data.fluxsect, t, output_data.total_flux, '*--');
            else
                % Fallback if dimensions don't match
                plot(t, output_data.fluxsect(1:length(t), :), t, output_data.total_flux(1:length(t)), '*--');
            end
            xlabel('Time [d]');
            ylabel('Sectional Flux [cm^3 m^{-2} d^{-1} sect^{-1}]');

            % Average settling velocity (subplot 2,4,8)
            subplot(2, 4, 8);
            % Match legacy exactly: plot(t_out, total_flux./total_mass/1e6)
            % Ensure dimensions match for plotting
            if length(t) == length(output_data.total_flux) && length(t) == length(output_data.total_mass)
                plot(t, output_data.total_flux ./ output_data.total_mass / 1e6);
            else
                % Fallback if dimensions don't match
                plot(t, output_data.total_flux(1:length(t)) ./ output_data.total_mass(1:length(t)) / 1e6);
            end
            xlabel('Time [d]');
            ylabel('Average v [m d^{-1}]');
        end

        function plotMassBalance(t, gains, losses)
            %PLOTMASSBALANCE Plot mass balance diagnostics
            figure(2);

            % Gains/Losses ratio (legacy Figure 2)
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
            %PLOTCONCENTRATIONEVOLUTION Plot concentration evolution
            figure(5);

            % Handle different dimension cases
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
            %PLOT3DSURFACES Plot 3D surface plots
            % fprintf('DEBUG: plot3DSurfaces called\n');
            n_sections = size(Y, 2);
            % fprintf('DEBUG: n_sections = %d, t length = %d, Y size = [%d,%d]\n', n_sections, length(t), size(Y,1), size(Y,2));

            % Safety checks for dimensions
            if length(t) ~= size(Y, 1) || n_sections == 0
                warning('Dimension mismatch in plot3DSurfaces: t=%d, Y=%dx%d', length(t), size(Y, 1), n_sections);
                return;
            end

            t_mat = t(:, ones(1, n_sections));
            % Fix: use v_lower (volume grid) to match legacy Figure 6
            if length(output_data.v_lower) == n_sections
                % Create v_mat as a matrix where each row is the same (v_lower values)
                % This matches legacy: v_mat = p.v_lower'; v_mat = v_mat(ones(length(t_out),1),:);
                v_mat = output_data.v_lower';
                v_mat = v_mat(ones(length(t), 1), :);
            else
                warning('Dimension mismatch: v_lower length=%d, n_sections=%d', length(output_data.v_lower), n_sections);
                v_mat = ones(length(t), n_sections); % fallback
            end

            % Gains/Losses ratio surface (legacy Figure 6) in Figure 4
            % Force exact match with legacy by constraining data range and colormap
            if size(gains, 1) == length(t) && size(gains, 2) == n_sections && ...
                    size(losses, 1) == length(t) && size(losses, 2) == n_sections

                % Calculate ratio with proper handling of edge cases
                ratio = gains ./ max(losses, eps);  % Avoid division by zero

                % Handle negative values and infinities properly
                ratio(losses <= 0) = NaN;  % Set to NaN where losses are non-positive
                ratio(~isfinite(ratio)) = NaN;  % Handle Inf and NaN values

                % Set negative ratios to NaN (exclude them from plot)
                ratio(ratio < 0) = NaN;  % Exclude negative values instead of converting to positive

                figure(4);
                surf(log10(v_mat), t_mat, ratio);
                xlabel('Volume');
                ylabel('Time [d]');
                zlabel('Gains/Losses');
                % title('Gains/Losses');

                % Set z-axis to start from 0, only considering positive values
                valid_ratios = ratio(~isnan(ratio) & ratio > 0);
                if ~isempty(valid_ratios)
                    zlim([0, max(valid_ratios)]);
                end
            end
        end

        function exportData(output_data, filename)
            %EXPORTDATA Export output data to file
            % output_data = struct from spectraAndFluxes
            % filename = output filename

            if nargin < 2
                filename = sprintf('coagulation_output_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
            end

            save(filename, 'output_data');
            fprintf('Output data saved to: %s\n', filename);
        end
    end
end
