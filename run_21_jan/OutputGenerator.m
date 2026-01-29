classdef OutputGenerator < handle
    %OUTPUTGENERATOR Generates outputs and visualizations for coagulation simulation
    %
    % NEW-2025-12-21:
    % - Column-mode support for flux spectrum field naming (fluxspec_i)
    % - Column-aware plot3DSurfaces to avoid Ns*Nz mismatch
    % - Column-aware plotAll routing (so you don't get blank/odd panels)
    %
    % NEW-2026-01-20:
    % - Add consistent "bin weight" for column inventory + export diagnostics:
    %     export_weight = "ones"  (number-like)
    %     export_weight = "vbin"  (volume-like)
    % - Keeps all old fields; adds new weighted fields so nothing breaks.
    %
    % NEW-2026-01-23 (CRITICAL UNIT FIX):
    % - Support state units: Y can be per m^3 OR per cm^3.
    % - If cfg.state_units == "m3", convert N from #/m^3 -> #/cm^3 using 1e-6.
    % - Adds cm^2-based flux/inventory fields to avoid 1e4 confusion:
    %     bottom_fluxsect_*_cm3cm2d, column_total_*_cm3cm2, etc.
    %
    % NEW-2026-01-23 (ADRIAN BUDGET READY):
    % - Adds "best" fields in BOTH m^2 and cm^2:
    %     total_mass_best_cm3cm2, total_flux_best_cm3cm2d, etc.
    % - Keeps plot compatibility (fluxspec_i remains per m^2/day).

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

                % ======================================================
                % NEW-2026-01-23: STATE UNITS HANDLING
                % Default keeps old behavior (assume per cm^3).
                % If config.state_units == "m3", convert to per cm^3.
                % ======================================================
                state_units = "cm3";
                try
                    if isprop(config,'state_units') && ~isempty(config.state_units)
                        state_units = string(config.state_units);
                    end
                catch
                end

                N_to_cm3 = 1.0;
                if strcmpi(state_units,"m3") || strcmpi(state_units,"per_m3")
                    N_to_cm3 = 1e-6;  % #/m^3 -> #/cm^3
                end

                % Radii/diameters
                r_i    = grid.getFractalRadii();
                r_v    = grid.getConservedRadii();
                diam_i = grid.getImageDiameters(config);
                diam_v = grid.getVolumeDiameters();

                % Settling velocities: [cm/s] -> [cm/day]
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

                % ======================================================
                % NEW-2026-01-20: choose a bin-weight for diagnostics
                %   - "ones": treat state as number-like per bin
                %   - "vbin": treat state as volume-like using bin volume weights
                % ======================================================
                weight_mode = "ones";
                try
                    if isprop(config,'export_weight') && ~isempty(config.export_weight)
                        weight_mode = string(config.export_weight);
                    end
                catch
                end

                % Build weights (Ns x 1)
                wbin = ones(Ns,1);
                try
                    if strcmpi(weight_mode,"vbin")
                        wtmp = [];
                        try
                            if ismethod(grid,'getBinVolumes')
                                wtmp = grid.getBinVolumes();
                            end
                        catch
                        end
                        if isempty(wtmp)
                            try
                                wtmp = Disaggregation.getBinVolumes_cm3(grid);
                            catch
                                wtmp = [];
                            end
                        end
                        if ~isempty(wtmp)
                            wbin = wtmp(:);
                        else
                            % fallback
                            wbin = ones(Ns,1);
                            weight_mode = "ones";
                        end
                    end
                catch
                    wbin = ones(Ns,1);
                    weight_mode = "ones";
                end

                % ======================================================
                % Inventory + export arrays
                % NOTE: We now store BOTH cm^2 and m^2 forms.
                % ======================================================

                % OLD fields kept (names unchanged)
                column_inventory_by_section_cm3m2 = zeros(n_times, Ns);
                column_total_inventory_cm3m2      = zeros(n_times, 1);

                bottom_fluxsect_cm3m2d   = zeros(n_times, Ns);
                bottom_total_flux_cm3m2d = zeros(n_times, 1);

                % NEW: weighted versions (same naming as you already used)
                column_inventory_by_section_weighted = zeros(n_times, Ns);
                column_total_inventory_weighted      = zeros(n_times, 1);

                bottom_fluxsect_weighted   = zeros(n_times, Ns);
                bottom_total_flux_weighted = zeros(n_times, 1);

                % NEW-2026-01-23: cm^2-based outputs (no 1e4 factor)
                column_inventory_by_section_cm3cm2 = zeros(n_times, Ns);
                column_total_inventory_cm3cm2      = zeros(n_times, 1);

                column_inventory_by_section_weighted_cm3cm2 = zeros(n_times, Ns);
                column_total_inventory_weighted_cm3cm2      = zeros(n_times, 1);

                bottom_fluxsect_cm3cm2d   = zeros(n_times, Ns);
                bottom_total_flux_cm3cm2d = zeros(n_times, 1);

                bottom_fluxsect_weighted_cm3cm2d   = zeros(n_times, Ns);
                bottom_total_flux_weighted_cm3cm2d = zeros(n_times, 1);

                % snapshots
                N_bottom  = zeros(n_times, Ns);
                N_surface = zeros(n_times, Ns);

                % legacy-like bottom spectra
                nspec_v_bottom   = zeros(n_times, Ns);
                masspec_v_bottom = zeros(n_times, Ns);

                % NEW: per-bin bottom flux (weighted) for size fractions
                Fbin_bottom_weighted = zeros(n_times, Ns);

                for it = 1:n_times
                    Nlayer_raw = squeeze(Y3(it, :, :));   % Nz x Ns

                    % NEW-2026-01-23: convert to per cm^3 if needed
                    Nlayer = Nlayer_raw * N_to_cm3;

                    ysurf = Nlayer(1, :);                 % 1 x Ns
                    ybot  = Nlayer(end, :);               % 1 x Ns

                    N_surface(it, :) = ysurf;
                    N_bottom(it, :)  = ybot;

                    % --------------------------------------------------
                    % OLD (kept): "inventory" computed from sum(Nlayer)
                    % These were intended as cm^3/m^2-style diagnostics.
                    % --------------------------------------------------
                    column_inventory_by_section_cm3cm2(it, :) = (sum(Nlayer, 1) * dz_cm);           % per cm^2
                    column_total_inventory_cm3cm2(it) = sum(column_inventory_by_section_cm3cm2(it, :));

                    column_inventory_by_section_cm3m2(it, :) = column_inventory_by_section_cm3cm2(it, :) * 1e4;  % per m^2
                    column_total_inventory_cm3m2(it) = sum(column_inventory_by_section_cm3m2(it, :));

                    % --------------------------------------------------
                    % NEW: weighted inventory (consistent with export weight)
                    % --------------------------------------------------
                    column_inventory_by_section_weighted_cm3cm2(it, :) = ( (sum(Nlayer, 1) .* (wbin(:)')) * dz_cm ); % per cm^2
                    column_total_inventory_weighted_cm3cm2(it) = sum(column_inventory_by_section_weighted_cm3cm2(it, :));

                    column_inventory_by_section_weighted(it, :) = column_inventory_by_section_weighted_cm3cm2(it, :) * 1e4; % per m^2
                    column_total_inventory_weighted(it) = sum(column_inventory_by_section_weighted(it, :));

                    % --------------------------------------------------
                    % OLD (kept): bottom flux (unweighted)
                    % --------------------------------------------------
                    bottom_fluxsect_cm3cm2d(it, :)  = (ybot .* set_vel_cmday(:)');    % per cm^2/day
                    bottom_total_flux_cm3cm2d(it)   = sum(bottom_fluxsect_cm3cm2d(it, :));

                    bottom_fluxsect_cm3m2d(it, :)   = bottom_fluxsect_cm3cm2d(it, :) * 1e4; % per m^2/day
                    bottom_total_flux_cm3m2d(it)    = sum(bottom_fluxsect_cm3m2d(it, :));

                    % --------------------------------------------------
                    % NEW: bottom flux (weighted)
                    % --------------------------------------------------
                    Fbin_cm3cm2d = (ybot(:) .* set_vel_cmday(:) .* wbin(:))';  % 1 x Ns, per cm^2/day
                    Fbin_bottom_weighted(it, :) = Fbin_cm3cm2d * 1e4;          % keep legacy "per m^2/day" version too

                    bottom_fluxsect_weighted_cm3cm2d(it, :)  = Fbin_cm3cm2d;
                    bottom_total_flux_weighted_cm3cm2d(it)   = sum(Fbin_cm3cm2d);

                    bottom_fluxsect_weighted(it, :)          = Fbin_cm3cm2d * 1e4; % per m^2/day
                    bottom_total_flux_weighted(it)           = sum(bottom_fluxsect_weighted(it, :));

                    % ---- legacy spectra on bottom layer (uses Nlayer per cm^3 now) ----
                    nspec_v_bottom(it, :)   = ybot ./ (1.5 * grid.v_lower') ./ grid.dwidth';
                    masspec_v_bottom(it, :) = ybot ./ grid.dwidth';
                end

                % cumulative export (safe)  (MUST be outside loop)
                if numel(t) < 2
                    export_cum_cm3m2    = zeros(size(t(:)));
                    export_cum_weighted = zeros(size(t(:)));
                    export_cum_cm3cm2    = zeros(size(t(:)));
                    export_cum_weighted_cm3cm2 = zeros(size(t(:)));
                else
                    export_cum_cm3m2    = cumtrapz(t(:), bottom_total_flux_cm3m2d(:));
                    export_cum_weighted = cumtrapz(t(:), bottom_total_flux_weighted(:));

                    export_cum_cm3cm2    = cumtrapz(t(:), bottom_total_flux_cm3cm2d(:));
                    export_cum_weighted_cm3cm2 = cumtrapz(t(:), bottom_total_flux_weighted_cm3cm2d(:));
                end

                % image/volume conversion ratio (legacy)
                diaratio = (config.fr_dim/3) * diam_v_mat ./ diam_i_mat;
                nspec_i_bottom   = nspec_v_bottom   .* diaratio;
                masspec_i_bottom = masspec_v_bottom .* diaratio;

                % ======================================================
                % NEW-2025-12-21 + NEW-2026-01-20:
                % Define fluxspec_i_bottom consistent with weighted export.
                % Keep it as per m^2/day for plotting compatibility.
                % ======================================================
                fluxspec_i_bottom = Fbin_bottom_weighted;  % [cm^3 m^-2 d^-1 per bin]

                % Assemble output
                output_data = struct();

                output_data.Y3 = Y3;
                output_data.Nz = Nz;
                output_data.Ns = Ns;
                output_data.z  = config.getZ();
                output_data.dz = config.dz;

                % NEW-2026-01-23: state unit info
                output_data.state_units = char(state_units);
                output_data.N_to_cm3 = N_to_cm3;

                % Closure-safe core fields (OLD kept)
                output_data.column_inventory_by_section_cm3m2 = column_inventory_by_section_cm3m2;
                output_data.column_total_inventory_cm3m2      = column_total_inventory_cm3m2;

                output_data.bottom_fluxsect_cm3m2d   = bottom_fluxsect_cm3m2d;
                output_data.bottom_total_flux_cm3m2d = bottom_total_flux_cm3m2d;
                output_data.export_cum_cm3m2         = export_cum_cm3m2;

                % NEW weighted fields (OLD names kept)
                output_data.export_weight_mode = weight_mode;
                output_data.wbin = wbin(:);

                output_data.column_inventory_by_section_weighted = column_inventory_by_section_weighted;
                output_data.column_total_inventory_weighted      = column_total_inventory_weighted;

                output_data.bottom_fluxsect_weighted   = bottom_fluxsect_weighted;
                output_data.bottom_total_flux_weighted = bottom_total_flux_weighted;
                output_data.export_cum_weighted        = export_cum_weighted;

                % NEW-2026-01-23: cm^2 versions
                output_data.column_inventory_by_section_cm3cm2 = column_inventory_by_section_cm3cm2;
                output_data.column_total_inventory_cm3cm2      = column_total_inventory_cm3cm2;

                output_data.column_inventory_by_section_weighted_cm3cm2 = column_inventory_by_section_weighted_cm3cm2;
                output_data.column_total_inventory_weighted_cm3cm2      = column_total_inventory_weighted_cm3cm2;

                output_data.bottom_fluxsect_cm3cm2d   = bottom_fluxsect_cm3cm2d;
                output_data.bottom_total_flux_cm3cm2d = bottom_total_flux_cm3cm2d;
                output_data.export_cum_cm3cm2         = export_cum_cm3cm2;

                output_data.bottom_fluxsect_weighted_cm3cm2d   = bottom_fluxsect_weighted_cm3cm2d;
                output_data.bottom_total_flux_weighted_cm3cm2d = bottom_total_flux_weighted_cm3cm2d;
                output_data.export_cum_weighted_cm3cm2         = export_cum_weighted_cm3cm2;

                % Backward compatibility fields
                output_data.bottom_total_flux = bottom_total_flux_cm3m2d;
                output_data.total_flux        = bottom_total_flux_cm3m2d;
                output_data.total_mass        = column_total_inventory_cm3m2; % inventory

                % Also provide compatibility "best" totals (weighted preferred)
                output_data.bottom_total_flux_best = bottom_total_flux_weighted;
                output_data.total_flux_best        = bottom_total_flux_weighted;
                output_data.total_mass_best        = column_total_inventory_weighted;

                % NEW-2026-01-23: "best" totals in cm^2 space (for Adrian budget)
                output_data.bottom_total_flux_best_cm3cm2d = bottom_total_flux_weighted_cm3cm2d;
                output_data.total_flux_best_cm3cm2d        = bottom_total_flux_weighted_cm3cm2d;
                output_data.total_mass_best_cm3cm2         = column_total_inventory_weighted_cm3cm2;

                % snapshots
                output_data.N_bottom  = N_bottom;
                output_data.N_surface = N_surface;

                % spectra (bottom layer)
                output_data.nspec_v   = nspec_v_bottom;
                output_data.masspec_v = masspec_v_bottom;
                output_data.nspec_i   = nspec_i_bottom;
                output_data.masspec_i = masspec_i_bottom;

                % NEW-2025-12-21: provide the missing field used by plotSpectraAndFluxes
                output_data.fluxspec_i = fluxspec_i_bottom;

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

                    diag.weight_mode = weight_mode;
                    diag.wbin = wbin;
                    diag.M_weighted = column_total_inventory_weighted;
                    diag.F_weighted = bottom_total_flux_weighted;
                    diag.ExportCum_weighted = export_cum_weighted;

                    diag.state_units = state_units;
                    diag.N_to_cm3 = N_to_cm3;
                    diag.M_cm3cm2 = column_total_inventory_cm3cm2;
                    diag.F_cm3cm2d = bottom_total_flux_cm3cm2d;

                    diag.M_best_cm3cm2 = output_data.total_mass_best_cm3cm2;
                    diag.F_best_cm3cm2d = output_data.total_flux_best_cm3cm2d;

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

        % --- everything else unchanged, but with column-aware routing ---
        function plotAll(t, Y, output_data, gains, losses, config, sectional_gains, sectional_losses, betas)

            % ==========================================================
            % NEW-2025-12-21: Column-aware plot routing
            % ==========================================================
            if isfield(output_data,'Nz') && isfield(output_data,'Ns') && ...
               isprop(config,'use_column') && config.use_column

                OutputGenerator.plotSpectraAndFluxes(t, Y, output_data);

                % NOTE: legacy "Total System Mass Balance" uses gains/losses
                % that are not strictly comparable in column mode (units).
                % Keep it, but it can look extreme / not interpretable.
                OutputGenerator.plotMassBalance(t, gains, losses);

                OutputGenerator.plotCoagVsSettRatio(t, losses);

                % Column-aware 3D plots (bottom layer only)
                OutputGenerator.plot3DSurfaces(t, Y, sectional_gains, sectional_losses, output_data);
                return;
            end

            % OLD (kept): slab routing
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

            % ==========================================================
            % NEW-2025-12-21: compatibility shim for renamed/missing fields
            % ==========================================================
            try
                if ~isfield(output_data,'fluxspec_i')
                    if isfield(output_data,'fluxspec_i_mat')
                        output_data.fluxspec_i = output_data.fluxspec_i_mat;
                    elseif isfield(output_data,'fluxspec')
                        output_data.fluxspec_i = output_data.fluxspec;
                    elseif isfield(output_data,'fluxspec_mat')
                        output_data.fluxspec_i = output_data.fluxspec_mat;
                    end
                end
            catch
            end

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

            % NEW-2025-12-21: force-create fluxspec_i if missing (debug-friendly)
            if ~isfield(output_data,'fluxspec_i')
                if isfield(output_data,'fluxspec_i_mat')
                    output_data.fluxspec_i = output_data.fluxspec_i_mat;
                elseif isfield(output_data,'fluxspec')
                    output_data.fluxspec_i = output_data.fluxspec;
                elseif isfield(output_data,'fluxspec_mat')
                    output_data.fluxspec_i = output_data.fluxspec_mat;
                else
                    disp('DEBUG: output_data fields are:');
                    disp(fieldnames(output_data));
                    error('Missing flux spectrum field. Add mapping based on fieldnames(output_data).');
                end
            end

            plot(output_data.diam_i_mat(:,2:end)', output_data.fluxspec_i(:,2:end)');
            xlabel('Particle image diameter [cm]');
            ylabel('Volume flux spectra [cm^3 m^{-2} d^{-1}]');
            axis tight;

            subplot(2, 4, 3);
            semilogy(t, Y, t, output_data.total_mass, '*--');
            xlabel('Time [d]');
            ylabel('Sectional concentration [vol/vol/sect]');

            subplot(2, 4, 7);
            if isfield(output_data,'fluxsect')
                if length(t) == size(output_data.fluxsect, 1) && length(t) == length(output_data.total_flux)
                    plot(t, output_data.fluxsect, t, output_data.total_flux, '*--');
                else
                    plot(t, output_data.fluxsect(1:length(t), :), t, output_data.total_flux(1:length(t)), '*--');
                end
                xlabel('Time [d]');
                ylabel('Sectional Flux [cm^3 m^{-2} d^{-1} sect^{-1}]');
            else
                % Column mode might not store fluxsect; keep subplot but avoid errors
                plot(0,0); axis([0 1 0 1]);
                title('fluxsect not stored (column)');
            end

            subplot(2, 4, 8);
            if length(t) == length(output_data.total_flux) && length(t) == length(output_data.total_mass)
                plot(t, output_data.total_flux ./ max(output_data.total_mass, eps) / 1e6);
            else
                plot(t, output_data.total_flux(1:length(t)) ./ max(output_data.total_mass(1:length(t)), eps) / 1e6);
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

        function plot3DSurfaces(t, Y, gains, losses, output_data)

            % ==========================================================
            % NEW-2025-12-21: Column-aware 3D plots (bottom layer only)
            % This prevents the Ns*Nz mismatch that caused blank/odd plots.
            % ==========================================================
            if isfield(output_data,'Nz') && isfield(output_data,'Ns') && isfield(output_data,'Y3')
                Ns = output_data.Ns;

                if length(t) ~= size(output_data.Y3, 1)
                    warning('Dimension mismatch in plot3DSurfaces (column): t=%d, Y3=%dx%dx%d', ...
                        length(t), size(output_data.Y3,1), size(output_data.Y3,2), size(output_data.Y3,3));
                    return;
                end

                % Build matrices
                t_mat = t(:, ones(1, Ns));
                v = output_data.v_lower(:);  % Ns x 1
                v_mat = log10(v(:)') .* ones(length(t),1);

                % Bottom layer concentration by section
                Nbot = squeeze(output_data.Y3(:, end, :)); % [time x Ns]

                % Plot bottom-layer concentration surface
                figure(4);
                clf;
                surf(v_mat, t_mat, Nbot, 'EdgeColor','none');
                xlabel('log10(Volume)');
                ylabel('Time [d]');
                zlabel('Bottom conc');
                title('Bottom-layer concentration surface (column)');
                view(35,20);

                % If gains/losses are provided as [time x Ns], also plot G/L surface
                if size(gains,1) == length(t) && size(gains,2) == Ns && ...
                   size(losses,1) == length(t) && size(losses,2) == Ns

                    ratio = gains ./ max(losses, eps);
                    ratio(losses <= 0) = NaN;
                    ratio(~isfinite(ratio)) = NaN;

                    figure(6);
                    clf;
                    surf(v_mat, t_mat, ratio, 'EdgeColor','none');
                    xlabel('log10(Volume)');
                    ylabel('Time [d]');
                    zlabel('Gains/Losses');
                    title('Bottom-layer Gains/Losses (column)');
                    view(35,20);
                end

                return;
            end

            % ==========================================================
            % ORIGINAL 0-D slab code (unchanged below)
            % ==========================================================
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