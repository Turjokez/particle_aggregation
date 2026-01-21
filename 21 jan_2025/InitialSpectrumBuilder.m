classdef InitialSpectrumBuilder < handle
    %INITIALSPECTRUMBUILDER Builds initial particle size spectra

    methods (Static)

        function spec_init = initialSpectrum(config, grid)
            %INITIALSPECTRUM Build initial spectrum with equal volume per section
            % Returns initial concentration vector for all sections
            %
            % NEW-2025-12-11:
            % If config.use_column = true, return a flattened column state
            % vector of size [Ns*Nz x 1], with a simple vertical profile.

            % ---------------------------
            % Original 0-D spectrum logic
            % ---------------------------
            spec0 = grid.av_vol(1) * ones(config.n_sections, 1);

            tfactor = 10.^(0:config.n_sections-1)';
            spec0 = max(spec0 ./ tfactor, 1.0e-30);

            spec0 = spec0 * config.num_1;

            % ----------------------------------------------------------
            % NEW-2025-12-11: Column option
            % ----------------------------------------------------------
            if isprop(config, 'use_column') && config.use_column
                Nz = config.getNumLayers();
                Ns = config.n_sections;

                % NEW-2025-12-11:
                % default vertical profile choice:
                %   - surface_only: all particles start in layer 1
                profile_type = 'surface_only';  % NEW-2025-12-11: default

                % Optional: allow config to define a profile type
                if isprop(config, 'init_profile') && ~isempty(config.init_profile)
                    profile_type = config.init_profile;  % NEW-2025-12-11
                end

                % Build Nz weights (sum of weights = 1)
                w_z = InitialSpectrumBuilder.columnProfile(config, profile_type);

                % Build Ns x Nz initial matrix: N(:,k) = spec0 * w_z(k)
                N0 = zeros(Ns, Nz);
                for k = 1:Nz
                    N0(:, k) = spec0 * w_z(k);
                end

                spec_init = N0(:);   % flatten to [Ns*Nz x 1]
                return;
            end

            % ---------------------------
            % Return 0-D default
            % ---------------------------
            spec_init = spec0;
        end


        % ==============================================================
        % NEW-2025-12-11: helper for vertical initial profile
        % ==============================================================
        function w_z = columnProfile(config, profile_type, varargin)
            %COLUMNPROFILE Build vertical weights for initial condition
            %
            % Returns w_z [Nz x 1] with sum(w_z)=1
            %
            % profile_type options:
            %   'surface_only' : all mass in layer 1
            %   'uniform'      : equal mass in all layers
            %   'exp'          : exponential decay with depth scale z0 (m)
            %
            % NEW-2025-12-12:
            %   'top_only'     : alias of surface_only 
            %
            % Optional args for 'exp':
            %   z0 (m): e-folding depth, default 30 m

            Nz = config.getNumLayers();
            z  = config.getZ();  % cell centers (m)

            if nargin < 2 || isempty(profile_type)
                profile_type = 'surface_only';
            end

            % ==========================================================
            % NEW-2025-12-12: alias support (do not break older scripts)
            % ==========================================================
            if ischar(profile_type) || isstring(profile_type)
                if strcmpi(profile_type, 'top_only')
                    profile_type = 'surface_only';
                end
            end

            switch lower(profile_type)

                case 'surface_only'
                    w_z = zeros(Nz, 1);
                    w_z(1) = 1.0;

                case 'uniform'
                    w_z = ones(Nz, 1) / Nz;

                case 'exp'
                    % NEW-2025-12-11: exponential decay with depth
                    if nargin >= 3
                        z0 = varargin{1};
                    elseif isprop(config, 'init_z0') && ~isempty(config.init_z0)
                        z0 = config.init_z0;          % NEW-2025-12-11
                    else
                        z0 = 30;                      % NEW-2025-12-11 default (m)
                    end

                    w_z = exp(-z(:) / z0);
                    if all(w_z == 0)
                        w_z = zeros(Nz, 1);
                        w_z(1) = 1.0;
                    else
                        w_z = w_z / sum(w_z);
                    end

                otherwise
                    error('Unknown column init profile_type: %s', profile_type);
            end
        end


        function spec_init = steadyStateSpectrum(config, grid, linear_ops)
            %STEADYSTATESPECTRUM Build steady-state spectrum (future implementation)
            % This would use fsolve to find steady state
            % For now, returns the simple initial spectrum

            warning('Steady state spectrum not yet implemented, using simple initial spectrum');
            spec_init = InitialSpectrumBuilder.initialSpectrum(config, grid);
        end


        function spec_init = customSpectrum(config, grid, spectrum_type, varargin)
            %CUSTOMSPECTRUM Build custom spectrum based on type
            % spectrum_type: 'uniform', 'exponential', 'power_law'

            switch spectrum_type
                case 'uniform'
                    spec_init = InitialSpectrumBuilder.uniformSpectrum(config, grid, varargin{:});
                case 'exponential'
                    spec_init = InitialSpectrumBuilder.exponentialSpectrum(config, grid, varargin{:});
                case 'power_law'
                    spec_init = InitialSpectrumBuilder.powerLawSpectrum(config, grid, varargin{:});
                otherwise
                    error('Unknown spectrum type: %s', spectrum_type);
            end
        end


        function spec_init = uniformSpectrum(config, grid, concentration)
            %UNIFORMSPECTRUM Build uniform spectrum across all sections
            if nargin < 3
                concentration = config.num_1;
            end

            spec_init = concentration * ones(config.n_sections, 1);
        end


        function spec_init = exponentialSpectrum(config, grid, base_concentration, decay_factor)
            %EXPONENTIALSPECTRUM Build exponentially decaying spectrum
            if nargin < 3
                base_concentration = config.num_1;
            end
            if nargin < 4
                decay_factor = 0.1;
            end

            spec_init = base_concentration * exp(-decay_factor * (0:config.n_sections-1)');
        end


        function spec_init = powerLawSpectrum(config, grid, base_concentration, power)
            %POWERLAWSPECTRUM Build power law spectrum
            if nargin < 3
                base_concentration = config.num_1;
            end
            if nargin < 4
                power = -2.0;
            end

            spec_init = base_concentration * (grid.av_vol / grid.av_vol(1)).^power;
        end
    end
end