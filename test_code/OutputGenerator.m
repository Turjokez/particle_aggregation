classdef OutputGenerator < handle
    %OUTPUTGENERATOR Generates outputs and visualizations for coagulation simulation

    methods (Static)
        function output_data = spectraAndFluxes(t, Y, gridobj, config)
            %SPECTRAANDFLUXES Compute particle spectra and fluxes
            % t = time vector
            % Y = concentration matrix
            % grid = DerivedGrid object
            % config = SimulationConfig object
            % Returns: struct with computed data

            n_times    = length(t);
            n_sections = length(gridobj.v_lower);

            % Initialize output arrays
            nspec_v   = zeros(n_times, n_sections);
            masspec_v = nspec_v;
            fluxsect  = nspec_v;
            fluxspec  = nspec_v;

            % Get radii and diameters
            r_i   = gridobj.getFractalRadii();
            r_v   = gridobj.getConservedRadii();
            diam_i = gridobj.getImageDiameters(config);
            diam_v = gridobj.getVolumeDiameters();

            % Calculate settling velocities
            set_vel = SettlingVelocityService.velocity(r_i, r_v, gridobj.setcon);
            set_vel = set_vel / 100 * config.day_to_sec;

            % Create time matrices for vectorized operations
            % diam_i and diam_v should be row vectors to match legacy
            diam_i     = diam_i';
            diam_v     = diam_v';
            diam_i_mat = diam_i(ones(n_times, 1), :);
            diam_v_mat = diam_v(ones(n_times, 1), :);

            % Compute spectra and fluxes for each time
            for jindx = 1:n_times
                yout = Y(jindx, :);

                nspec_v(jindx, :)   = yout ./ (1.5 * gridobj.v_lower') ./ gridobj.dwidth';
                masspec_v(jindx, :) = yout ./ gridobj.dwidth';
                fluxsect(jindx, :)  = yout .* set_vel' * 1e6;
                fluxspec(jindx, :)  = masspec_v(jindx, :) .* set_vel' * 1e6;
            end

            % Total quantities
            total_flux = sum(fluxsect, 2);
            total_mass = sum(Y, 2);

            % Diameter ratio and image-based spectra
            diaratio  = (config.fr_dim/3) * diam_v_mat ./ diam_i_mat;
            nspec_i   = nspec_v .* diaratio;
            masspec_i = masspec_v .* diaratio;
            fluxspec_i = fluxspec .* diaratio;

            % Assemble output data
            output_data                 = struct();
            output_data.nspec_v         = nspec_v;
            output_data.masspec_v       = masspec_v;
            output_data.fluxsect        = fluxsect;
            output_data.fluxspec        = fluxspec;
            output_data.total_flux      = total_flux;
            output_data.total_mass      = total_mass;
            output_data.diam_i          = diam_i;
            output_data.diam_v          = diam_v;
            output_data.diam_i_mat      = diam_i_mat;
            output_data.diam_v_mat      = diam_v_mat;
            output_data.nspec_i         = nspec_i;
            output_data.masspec_i       = masspec_i;
            output_data.fluxspec_i      = fluxspec_i;
            output_data.set_vel         = set_vel;
            output_data.v_lower         = gridobj.v_lower;
            output_data.dwidth          = gridobj.dwidth;
            output_data.diaratio        = diaratio;
            output_data.fr_dim          = config.fr_dim;

            % Diagnostics snapshot for cross-checking legacy vs OOP
            try
                diag = struct();
                diag.t        = t;
                diag.Y        = Y;
                diag.r_i      = r_i;
                diag.r_v      = r_v;
                diag.diam_i   = diam_i;
                diag.diam_v   = diam_v;
                diag.set_vel  = set_vel;
                diag.v_lower  = gridobj.v_lower;
                diag.dwidth   = gridobj.dwidth;
                diag.nspec_v  = nspec_v;
                diag.masspec_v= masspec_v;
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

            % Figure 1: Spectra and fluxes (match legacy axes/layout)
            OutputGenerator.plotSpectraAndFluxes(t, Y, output_data);

            % Figure 2: Mass balance (single ratio panel only, like legacy)
            OutputGenerator.plotMassBalance(t, gains, losses);

            % Figure 3: (Coag Losses)/(Settling Losses) vs time (legacy)
            OutputGenerator.plotCoagVsSettRatio(t, losses);

            % Figure 4: Legacy Figure 6 → sectional (Gains/Losses) surface
            OutputGenerator.plot3DSurfaces(t, Y, sectional_gains, sectional_losses, output_data);
        end

        function plotCoagVsSettRatio(t, losses)
            % Legacy Figure 3: (Coag Losses)/(Settling Losses)
            figure(3); clf;
        
            % pick whichever settling-loss field exists
            if ~isfield(losses,'coag') || ~(isfield(losses,'sett') || isfield(losses,'settl'))
                warning('OutputGenerator:plotCoagVsSettRatio', ...
                    'Missing losses.coag or settling field (sett/settl). Skipping plot.');
                return;
            end
        
            if isfield(losses,'sett')
                den = losses.sett;
            else
                den = losses.settl;
            end
        
            % vectorize + guardrails
            t    = t(:);
            num  = losses.coag(:);
            den  = den(:);
        
            % align lengths if needed
            n = min([numel(t) numel(num) numel(den)]);
            t   = t(1:n);
            num = num(1:n);
            den = max(den(1:n), eps);     % avoid divide-by-zero/tiny
        
            ratio = num ./ den;
            ratio(~isfinite(ratio)) = NaN;
        
            plot(t, ratio, 'LineWidth', 1.5);
            grid on;
            set(gca, 'FontName', 'Helvetica', 'FontSize', 14);
            xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 14);
            ylabel('(Coag Losses)/(Settling Losses)', 'FontName', 'Helvetica', 'FontSize', 14);
            title('Coagulation vs Settling Losses');
        end

        function plotSpectraAndFluxes(t, Y, output_data)
            %PLOTSPECTRAANDFLUXES Plot particle spectra and fluxes
            % Match legacy layout exactly: mixed 2x2 and 2x4 subplot layout
            figure(1);

            % --- GUARDRAIL: temporarily silence log-axes warning
            warnID = 'MATLAB:Axes:NegativeDataInLogAxes';
            prevW  = warning('query', warnID);
            warning('off', warnID);
            c = onCleanup(@() warning(prevW.state, warnID));

            % Number spectrum (subplot 2,2,1 - overlaid on 2x4 layout)
            subplot(2, 2, 1);
            ns_init  = output_data.nspec_i(1, :);
            ns_final = output_data.nspec_i(end, :);

            % --- GUARDRAIL: only positive/finite values for loglog
            [x1,y1] = OutputGenerator.pair_pos(output_data.diam_i(:), ns_init(:));
            [x2,y2] = OutputGenerator.pair_pos(output_data.diam_i(:), ns_final(:));

            loglog(x1, y1, 'b', x2, y2, 'r');
            xlabel('Particle diameter [cm]');
            ylabel('Number spectrum [# cm^{-4}]');
            axis tight;

            % Flux spectra (subplot 2,2,3 - overlaid on 2x4 layout) — linear axes
            subplot(2, 2, 3);
            fluxspec_i_plot = max(output_data.fluxspec_i, 0);  % clip negatives for plotting
            plot(output_data.diam_i_mat(:,2:end)', fluxspec_i_plot(:,2:end)');
            xlabel('Particle image diameter [cm]');
            ylabel('Volume flux spectra [cm^2 m^{-2} d^{-1}]');
            axis tight;

            % Sectional concentration vs time (subplot 2,4,3)
            subplot(2, 4, 3);
            Yplot = max(Y, eps);                      % --- GUARDRAIL
            TM    = max(output_data.total_mass, eps); % --- GUARDRAIL
            semilogy(t, Yplot, t, TM, '*--');
            xlabel('Time [d]');
            ylabel('Sectional concentration [vol/vol/sect]');

            % Sectional flux over time with total (subplot 2,4,7)
            subplot(2, 4, 7);
            if length(t) == size(output_data.fluxsect, 1) && length(t) == length(output_data.total_flux)
                plot(t, output_data.fluxsect, t, output_data.total_flux, '*--');
            else
                plot(t, output_data.fluxsect(1:length(t), :), t, output_data.total_flux(1:length(t)), '*--');
            end
            xlabel('Time [d]');
            ylabel('Sectional Flux [cm^3 m^{-2} d^{-1} sect^{-1}]');

            % Average settling velocity (subplot 2,4,8)
            subplot(2, 4, 8);
            if length(t) == length(output_data.total_flux) && length(t) == length(output_data.total_mass)
                vbar = output_data.total_flux ./ max(output_data.total_mass, eps) / 1e6; % --- GUARDRAIL
                plot(t, vbar);
            else
                vbar = output_data.total_flux(1:length(t)) ./ max(output_data.total_mass(1:length(t)), eps) / 1e6; % --- GUARDRAIL
                plot(t, vbar);
            end
            xlabel('Time [d]');
            ylabel('Average v [m d^{-1}]');
        end

        function plotMassBalance(t, gains, losses)
            %PLOTMASSBALANCE Plot mass balance diagnostics
            figure(2); clf;

            if isfield(gains, 'growth') && isfield(losses, 'sett') && isfield(losses, 'coag')
                denominator = losses.sett + losses.coag;
                ratio = gains.growth ./ max(denominator, eps);   % --- GUARDRAIL
                ratio(~isfinite(ratio)) = NaN;                   % --- GUARDRAIL
                plot(t, ratio);
                xlabel('Time [d]');
                ylabel('Gains/Losses');
                title('Total System Mass Balance');
            end
        end

        function plotConcentrationEvolution(t, Y, total_mass)
            %PLOTCONCENTRATIONEVOLUTION Plot concentration evolution
            figure(5);

            if size(Y, 1) == length(t) && size(Y, 2) > 0
                semilogy(t, max(Y, eps));  % --- GUARDRAIL
                hold on;
                if length(total_mass) == length(t)
                    semilogy(t, max(total_mass, eps), '*--', 'LineWidth', 2); % --- GUARDRAIL
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
            %PLOT3DSURFACES Plot 3D surface plots (legacy Figure 6 styling)
            n_sections = size(Y, 2);

            if length(t) ~= size(Y, 1) || n_sections == 0
                warning('Dimension mismatch in plot3DSurfaces: t=%d, Y=%dx%d', length(t), size(Y, 1), n_sections);
                return;
            end

            t_mat = t(:, ones(1, n_sections));

            % Use v_lower (volume grid) to match legacy
            if length(output_data.v_lower) == n_sections
                v_mat = output_data.v_lower';
                v_mat = v_mat(ones(length(t), 1), :);
            else
                warning('Dimension mismatch: v_lower length=%d, n_sections=%d', length(output_data.v_lower), n_sections);
                v_mat = ones(length(t), n_sections); % fallback
            end

            % Gains/Losses ratio surface (Figure 4)
            if size(gains,1)==length(t) && size(gains,2)==n_sections && ...
               size(losses,1)==length(t) && size(losses,2)==n_sections

                % --- PATCH: show (coag losses)/(settling losses) per section (like Fig. 3)
                if isstruct(losses) && isfield(losses,'coag') && isfield(losses,'sett')
                    denom = max(losses.sett, eps);
                    ratio = losses.coag ./ denom;
                else
                    % fallback if 'losses' isn't split into fields
                    ratio = gains ./ max(losses, eps);
                end
                ratio(~isfinite(ratio) | ratio <= 0) = NaN;   % keep guardrail

                figure(4);
                surf(log10(v_mat), t_mat, ratio);
                xlabel('Volume');
                ylabel('Time [d]');
                zlabel('Coag/Sett');

                valid_ratios = ratio(~isnan(ratio) & ratio > 0);
                if ~isempty(valid_ratios)
                    zlim([0, max(valid_ratios)]);
                end
            end
        end

        function exportData(output_data, filename)
            %EXPORTDATA Export output data to file
            if nargin < 2
                filename = sprintf('coagulation_output_%s.mat', datestr(now, 'yyyymmdd_HHMMSS'));
            end

            save(filename, 'output_data');
            fprintf('Output data saved to: %s\n', filename);
        end

        %% Flux vs Time plotting helper ===
        function [t, flux] = simpleFluxTimeSeries(result, gridobj)
            % SIMPLEFLUXTIMESERIES  Compute total flux vs. time
            t = result.time(:);                 % Time vector
            Y = result.concentrations;          % Concentration matrix (time × sections)

            % Settling velocity or proxy weighting
            if isfield(result, 'w_sec')
                w = result.w_sec(:)';            % Settling speed per section
            else
                n   = size(Y, 2);
                r_i = gridobj.getFractalRadii();
                r_v = gridobj.getConservedRadii();
                set_vel = SettlingVelocityService.velocity(r_i, r_v, gridobj.setcon);
                set_vel = set_vel(:)' / 100 * 86400;  % cm/s → m/day
                w = set_vel / max(set_vel);           % normalized velocity profile
            end

            flux = sum(Y .* w, 2);

            figure('Color','w');
            plot(t, flux, 'b', 'LineWidth', 2);
            xlabel('Time (days)');
            ylabel('Total flux (a.u.)');
            title('Flux vs Time');
            grid on;
        end

        %% === Step 5: Flux vs Size vs Time Surface ===
        function plotFluxSurface(result, gridobj)
            % PLOTFLUXSURFACE  3-D surface of flux vs. size vs. time
            t = result.time(:);
            Y = result.concentrations;

            r_i = gridobj.getFractalRadii();
            r_v = gridobj.getConservedRadii();
            set_vel = SettlingVelocityService.velocity(r_i, r_v, gridobj.setcon);
            set_vel = set_vel(:)' / 100 * 86400;   % cm/s → m d⁻¹

            flux_section = Y .* set_vel;

            figure('Color','w');
            flux_section = max(flux_section, 0);   % --- GUARDRAIL
            
            D = max(2*r_i, realmin);              % --- GUARDRAIL (strictly positive)
            [T, D] = meshgrid(t, D);
            surf(T', log10(D'), flux_section, 'EdgeColor','none');
            xlabel('Time [days]');
            ylabel('log_{10} Diameter [cm]');
            zlabel('Flux [a.u.]');
            title('Flux vs Size vs Time');
            colorbar;
            view(45,30);
            grid on;
        end   % closes plotFluxSurface
               
        %% === Step 7: Gains / Losses Surface Diagnostic ===
        function plotGainsLossesSurface(result, output_data)
            % PLOTGAINSLOSSESSURFACE  Visualize Gains/Losses ratio in 3D
            fprintf('Generating Gains/Losses diagnostic surface...\n');
            
            if ~isfield(result, 'diagnostics') || ...
               ~isfield(result.diagnostics, 'gains') || ...
               ~isfield(result.diagnostics, 'losses')
                warning('No gain/loss diagnostics found in result struct.');
                return;
            end
            
            gains  = result.diagnostics.gains;
            losses = result.diagnostics.losses;
            t      = result.time(:);
            v      = output_data.v_lower;
            
            if size(gains,1) ~= length(t)
                warning('Dimension mismatch: gains(%d) vs time(%d)', size(gains,1), length(t));
                return;
            end
            
            ratio = gains ./ max(losses, eps);      % --- GUARDRAIL
            ratio(ratio <= 0) = NaN;                % --- GUARDRAIL
            ratio(~isfinite(ratio)) = NaN;          % --- GUARDRAIL
            
            [T, V] = meshgrid(t, log10(v));
            
            figure('Color','w');
            surf(T', V', ratio, 'EdgeColor','none');
            xlabel('Time [d]');
            ylabel('log_{10} Volume [cm^3]');
            zlabel('Gains / Losses');
            title('Step 7 — Coagulation vs Settling Surface');
            colorbar;
            view(45,30);
            grid on;
        end
    end

    % --- tiny private helpers (so you don't need a separate utility file)
    methods (Static, Access = private)
        function y = pos(y, tiny)
            if nargin < 2 || isempty(tiny), tiny = 1e-30; end
            y = double(y);
            y(~isfinite(y) | y <= tiny) = NaN;
        end
        function [x2,y2] = pair_pos(x, y, tiny)
            if nargin < 3, tiny = 1e-30; end
            x = x(:); y = y(:);
            ok = isfinite(x) & isfinite(y) & (y > tiny);
            if any(ok), x2 = x(ok); y2 = y(ok);
            else,       x2 = NaN;   y2 = NaN;
            end
        end
    end
end