classdef DerivedGrid < handle
    %DERIVEDGRID Precomputed constants and grids derived from configuration
    %
    % NEW-2025-12-11:
    % Adds a simple vertical grid for 1-D column experiments.
    % Keeps original 0-D slab behavior when config.use_column = false.

    properties
        % ---------------------------
        % Derived constants
        % ---------------------------
        dvisc;          % Dynamic viscosity [g cm^{-1} s^{-1}]
        del_rho;        % Density difference parameter
        conBr;          % Brownian constant
        a0;             % Unit particle radius [cm]
        v0;             % Unit particle volume [cm^3]

        % ---------------------------
        % Grid arrays (size bins)
        % ---------------------------
        v_lower;        % Lower volume bounds for sections
        v_upper;        % Upper volume bounds for sections
        av_vol;         % Average volume per section
        dcomb;          % Combined diameter
        dwidth;         % Section width

        % ---------------------------
        % Fractal parameters
        % ---------------------------
        amfrac;         % Fractal radius multiplier
        bmfrac;         % Fractal dimension exponent

        % ---------------------------
        % Settling parameters
        % ---------------------------
        setcon;         % Settling velocity constant

        % ==========================================================
        % NEW-2025-12-11: VERTICAL GRID FOR 1-D COLUMN
        % ==========================================================
        z_edges;        % layer interfaces [m], length Nz+1
        z_centers;      % layer centers [m], length Nz
        Nz;             % number of vertical layers
    end

    methods
        function obj = DerivedGrid(config)
            %DERIVEDGRID Constructor - compute derived values from config
            obj.computeDerivedValues(config);
        end

        function computeDerivedValues(obj, config)
            %COMPUTEDERIVEDVALUES Compute all derived constants and grids

            % ---------------------------
            % Dynamic viscosity
            % ---------------------------
            obj.dvisc = config.kvisc * config.rho_fl;

            % ---------------------------
            % Density difference parameter
            % ---------------------------
            obj.del_rho = (4.5*2.48) * config.kvisc * config.rho_fl / config.g * (config.d0/2)^(-0.83);

            % ---------------------------
            % Brownian constant
            % ---------------------------
            obj.conBr = 2.0/3.0 * config.k * config.temp / obj.dvisc;

            % ---------------------------
            % Unit particle properties
            % ---------------------------
            obj.a0 = config.d0 / 2;
            obj.v0 = (pi/6) * config.d0^3;

            % ---------------------------
            % Volume grid (size sections)
            % ---------------------------
            obj.v_lower = obj.v0 * 2.^(0:config.n_sections-1)';
            obj.v_upper = 2.0 * obj.v_lower;
            obj.av_vol  = 1.5 * obj.v_lower;

            % ---------------------------
            % Diameter grid
            % ---------------------------
            obj.dcomb  = (obj.v_lower * 6/pi).^(1.0/3.0);
            obj.dwidth = (2^(1.0/3.0) - 1) * obj.dcomb;

            % ---------------------------
            % Fractal parameters
            % ---------------------------
            amfrac_temp = (4.0/3.0*pi)^(-1.0/config.fr_dim) * obj.a0^(1.0-3.0/config.fr_dim);
            obj.amfrac  = amfrac_temp * sqrt(0.6);
            obj.bmfrac  = 1.0/config.fr_dim;

            % ---------------------------
            % Settling constant
            % ---------------------------
            obj.setcon = (2.0/9.0) * obj.del_rho / config.rho_fl * config.g / config.kvisc;

            % ======================================================
            % NEW-2025-12-11: VERTICAL GRID (SLAB OR 1-D COLUMN)
            % ======================================================
            %
            % OLD SLAB INTERPRETATION (kept):
            %   treat the domain as one layer of thickness dz
            %
            % NEW COLUMN MODE:
            %   Nz = floor(z_max / dz)
            %   z_edges   = 0:dz:(Nz*dz)
            %   z_centers = midpoint of each layer
            %
            % NOTE:
            %   We do NOT force z_edges(end) to equal z_max exactly.
            %   It becomes Nz*dz. That is consistent with getNumLayers().
            %
            % (If you want exact z_max later, we can add a "last partial layer"
            %  option — but for now keep it simple and stable.)

            if isprop(config, 'use_column') && config.use_column
                z_max = config.z_max;  % [m]
                dz    = config.dz;     % [m]

                if dz <= 0
                    error('DerivedGrid: dz must be positive when use_column = true.');
                end
                if z_max <= 0
                    error('DerivedGrid: z_max must be positive when use_column = true.');
                end

                obj.Nz = floor(z_max / dz);
                if obj.Nz < 1
                    obj.Nz = 1;
                end

                % OLD (kept, but commented out): row vector
                % obj.z_edges   = (0:obj.Nz) * dz;  % [m]

                % NEW: column vector (consistent with slab case)
                obj.z_edges   = ((0:obj.Nz)' * dz);  % [m]

                obj.z_centers = (obj.z_edges(1:end-1) + obj.z_edges(2:end)) / 2;

            else
                % Default slab: one layer from 0 to dz
                obj.Nz = 1;
                obj.z_edges   = [0; config.dz];
                obj.z_centers = config.dz / 2;
            end
        end

        function r_i = getFractalRadii(obj)
            %GETFRACTALRADII Get fractal radii for all sections
            r_i = obj.amfrac * obj.av_vol.^obj.bmfrac;
        end

        function r_v = getConservedRadii(obj)
            %GETCONSERVEDRADII Get conserved volume radii for all sections
            r_v = (0.75/pi * obj.av_vol).^(1.0/3.0);
        end

        function d_i = getImageDiameters(obj, config)
            %GETIMAGEDIAMETERS Get image diameters for all sections
            r_i = obj.getFractalRadii();
            d_i = 2.0 * config.r_to_rg * r_i;
        end

        function d_v = getVolumeDiameters(obj)
            %GETVOLUMEDIAMETERS Get volume diameters for all sections
            r_v = obj.getConservedRadii();
            d_v = 2.0 * r_v;
        end
    end
end
