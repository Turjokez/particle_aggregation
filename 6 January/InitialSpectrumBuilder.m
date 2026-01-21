% ==============================================================
% FILE: InitialSpectrumBuilder.m
% ==============================================================
% UPDATED (2026-01-14):
% - FIX: ensure ONLY ONE classdef per file (MATLAB requirement).
% - Builds NUMBER spectrum n0 first, then converts to BIOVOLUME if requested.
% - NEW SAFETY RESCALE: forces n0(1) == cfg.num_1 (preferred) or cfg.n1.
%   This prevents catastrophic IC magnitudes (e.g., 1e11 in bin 1).
% - Column mode: applies vertical weights w_z with sum(w_z)=1.
% - Adds optional cfg.debug_init_spectrum prints (default false).
% ==============================================================

classdef InitialSpectrumBuilder < handle
    %INITIALSPECTRUMBUILDER Builds initial particle size spectra

    methods (Static)

        function spec_init = initialSpectrum(config, grid)
            %INITIALSPECTRUM Build initial spectrum
            %
            % IMPORTANT (NEW-2026-01-13):
            %   Respect state convention:
            %     - state_is_biovolume = false: return NUMBER spectrum n0
            %     - state_is_biovolume = true : return BIOVOLUME spectrum v0 = n0 .* av_vol
            %
            % Column option:
            %   if use_column=true, returns flattened [Ns*Nz x 1] with vertical profile weights.

            Ns = config.n_sections;

            % ---------------------------
            % 0) Debug flag (optional)
            % ---------------------------
            dbg = false;
            try
                if isprop(config,'debug_init_spectrum') && ~isempty(config.debug_init_spectrum)
                    dbg = logical(config.debug_init_spectrum);
                end
            catch
                dbg = false;
            end

            % ---------------------------
            % 1) Build NUMBER spectrum shape (dimensionless)
            % ---------------------------
            % Default: decreasing by decades with bin index (your old tfactor logic)
            nshape = ones(Ns, 1);
            tfactor = 10.^(0:Ns-1)';     % 1,10,100,...
            nshape = nshape ./ tfactor;
            nshape = max(nshape, 1.0e-30);

            % ---------------------------
            % 2) Choose amplitude target for bin-1
            % ---------------------------
            amp = NaN;
            if isprop(config,'num_1') && ~isempty(config.num_1) && isfinite(config.num_1) && config.num_1 > 0
                amp = config.num_1;
            elseif isprop(config,'n1') && ~isempty(config.n1) && isfinite(config.n1) && config.n1 > 0
                amp = config.n1;
            else
                amp = 1e3; % safe fallback
            end

            % ---------------------------
            % 3) Apply amplitude in a SAFE way
            % ---------------------------
            % OLD/LEGACY (kept as comment):
            % n0 = nshape * amp;
            %
            % NEW SAFETY:
            % force n0(1) == amp exactly, regardless of shape definition.
            n0 = nshape;
            if isfinite(n0(1)) && n0(1) > 0
                n0 = n0 * (amp / n0(1));
            else
                n0 = nshape * amp; % fallback
            end

            n0 = max(n0, 0); % no negatives

            if dbg
                fprintf('[initSpectrum] Ns=%d, amp(bin1)=%.3e, n0(1)=%.3e, n0 med=%.3e, n0 max=%.3e\n', ...
                    Ns, amp, n0(1), median(n0(n0>0)), max(n0));
            end

            % ---------------------------
            % 4) Convert to state variable (n or v)
            % ---------------------------
            state_is_bv = false;
            if isprop(config,'state_is_biovolume') && ~isempty(config.state_is_biovolume)
                state_is_bv = logical(config.state_is_biovolume);
            end

            if state_is_bv
                av = grid.av_vol(:);
                assert(numel(av)==Ns, 'grid.av_vol length mismatch with n_sections');
                spec0 = n0 .* av;     % BIOVOLUME state
            else
                spec0 = n0;           % NUMBER state
            end

            % ----------------------------------------------------------
            % 5) Column option
            % ----------------------------------------------------------
            if isprop(config, 'use_column') && ~isempty(config.use_column) && config.use_column
                Nz = config.getNumLayers();

                profile_type = 'surface_only';
                if isprop(config, 'init_profile') && ~isempty(config.init_profile)
                    profile_type = config.init_profile;
                end

                w_z = InitialSpectrumBuilder.columnProfile(config, profile_type);

                % safety normalize (even if profile already does it)
                sw = sum(w_z);
                if ~isfinite(sw) || sw <= 0
                    w_z = zeros(Nz,1); w_z(1)=1;
                else
                    w_z = w_z / sw;
                end

                N0 = zeros(Ns, Nz);
                for k = 1:Nz
                    N0(:, k) = spec0 * w_z(k);
                end

                spec_init = N0(:);

                if dbg
                    fprintf('[initSpectrum] column Nz=%d, sum(w_z)=%.6f, sum(spec_init)=%.3e\n', ...
                        Nz, sum(w_z), sum(spec_init));
                end
                return;
            end

            % ---------------------------
            % Return 0-D default
            % ---------------------------
            spec_init = spec0;
        end


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
            % 'top_only' is alias of 'surface_only'

            Nz = config.getNumLayers();
            z  = config.getZ();  % cell centers (m)

            if nargin < 2 || isempty(profile_type)
                profile_type = 'surface_only';
            end

            if ischar(profile_type) || isstring(profile_type)
                if strcmpi(profile_type, 'top_only')
                    profile_type = 'surface_only';
                end
            end

            switch lower(string(profile_type))

                case "surface_only"
                    w_z = zeros(Nz, 1);
                    w_z(1) = 1.0;

                case "uniform"
                    w_z = ones(Nz, 1) / Nz;

                case "exp"
                    if nargin >= 3
                        z0 = varargin{1};
                    elseif isprop(config, 'init_z0') && ~isempty(config.init_z0)
                        z0 = config.init_z0;
                    else
                        z0 = 30;
                    end

                    w_z = exp(-z(:) / z0);
                    if all(w_z == 0) || ~any(isfinite(w_z))
                        w_z = zeros(Nz, 1);
                        w_z(1) = 1.0;
                    else
                        w_z = w_z / sum(w_z);
                    end

                otherwise
                    error('Unknown column init profile_type: %s', string(profile_type));
            end
        end


        function spec_init = steadyStateSpectrum(config, grid, linear_ops) %#ok<INUSD>
            warning('Steady state spectrum not yet implemented, using simple initial spectrum');
            spec_init = InitialSpectrumBuilder.initialSpectrum(config, grid);
        end


        function spec_init = customSpectrum(config, grid, spectrum_type, varargin)
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


        function spec_init = uniformSpectrum(config, grid, concentration) %#ok<INUSD>
            if nargin < 3
                concentration = config.num_1;
            end
            spec_init = concentration * ones(config.n_sections, 1);
        end


        function spec_init = exponentialSpectrum(config, grid, base_concentration, decay_factor) %#ok<INUSD>
            if nargin < 3
                base_concentration = config.num_1;
            end
            if nargin < 4
                decay_factor = 0.1;
            end
            spec_init = base_concentration * exp(-decay_factor * (0:config.n_sections-1)');
        end


        function spec_init = powerLawSpectrum(config, grid, base_concentration, power)
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