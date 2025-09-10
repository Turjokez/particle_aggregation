classdef DerivedGrid < handle
    %DERIVEDGRID Precomputed constants and grids derived from configuration
    
    properties
        % Derived constants
        dvisc;          % Dynamic viscosity [g cm^{-1} s^{-1}]
        del_rho;        % Density difference parameter
        conBr;          % Brownian constant
        a0;             % Unit particle radius [cm]
        v0;             % Unit particle volume [cm^3]
        
        % Grid arrays
        v_lower;        % Lower volume bounds for sections
        v_upper;        % Upper volume bounds for sections
        av_vol;         % Average volume per section
        dcomb;          % Combined diameter
        dwidth;         % Section width
        
        % Fractal parameters
        amfrac;         % Fractal radius multiplier
        bmfrac;         % Fractal dimension exponent
        
        % Settling parameters
        setcon;         % Settling velocity constant
    end
    
    methods
        function obj = DerivedGrid(config)
            %DERIVEDGRID Constructor - compute derived values from config
            obj.computeDerivedValues(config);
        end
        
        function computeDerivedValues(obj, config)
            %COMPUTEDERIVEDVALUES Compute all derived constants and grids
            
            % Dynamic viscosity
            obj.dvisc = config.kvisc * config.rho_fl;
            
            % Density difference parameter
            obj.del_rho = (4.5*2.48) * config.kvisc * config.rho_fl / config.g * (config.d0/2)^(-0.83);
            
            % Brownian constant
            obj.conBr = 2.0/3.0 * config.k * config.temp / obj.dvisc;
            
            % Unit particle properties
            obj.a0 = config.d0 / 2;
            obj.v0 = (pi/6) * config.d0^3;
            
            % Volume grid
            obj.v_lower = obj.v0 * 2.^(0:config.n_sections-1)';
            obj.v_upper = 2.0 * obj.v_lower;
            obj.av_vol = 1.5 * obj.v_lower;
            
            % Diameter grid
            obj.dcomb = (obj.v_lower * 6/pi).^(1.0/3.0);
            obj.dwidth = (2^(1.0/3.0) - 1) * obj.dcomb;
            
            % Fractal parameters
            amfrac_temp = (4.0/3.0*pi)^(-1.0/config.fr_dim) * obj.a0^(1.0-3.0/config.fr_dim);
            obj.amfrac = amfrac_temp * sqrt(0.6);
            obj.bmfrac = 1.0/config.fr_dim;
            
            % Settling constant
            obj.setcon = (2.0/9.0) * obj.del_rho / config.rho_fl * config.g / config.kvisc;
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
