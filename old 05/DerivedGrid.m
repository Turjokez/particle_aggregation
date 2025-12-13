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
            
            % --- Dynamic viscosity [g cm^-1 s^-1]
            obj.dvisc = config.kvisc * config.rho_fl;
            
            % --- Density difference parameter
            obj.del_rho = (4.5 * 2.48) * config.kvisc * config.rho_fl / ...
                          config.g * (config.d0/2)^(-0.83);
            
            % --- Brownian constant
            obj.conBr = (2.0/3.0) * config.k * config.temp / obj.dvisc;
            
            % --- Unit particle properties
            obj.a0 = config.d0 / 2;
            obj.v0 = (pi/6) * config.d0^3;
            
            % --- Volume grid
            obj.v_lower = obj.v0 * 2.^(0:config.n_sections-1)';
            obj.v_upper = 2.0 * obj.v_lower;
            obj.av_vol  = 1.5 * obj.v_lower;
            
            % ---- Debug: report smallest & largest diameters ---------------
            d_min = (obj.v_lower(1) * 6/pi)^(1/3);
            d_max = (obj.v_upper(end) * 6/pi)^(1/3);
            fprintf('DerivedGrid → diameter range: %.3e–%.3e cm\n', d_min, d_max);
            
            % --- Diameter grid
            obj.dcomb  = (obj.v_lower * 6/pi).^(1.0/3.0);
            obj.dwidth = (2^(1.0/3.0) - 1) * obj.dcomb;
            
            % --- Fractal parameters
            amfrac_temp = (4.0/3.0*pi)^(-1.0/config.fr_dim) * ...
                          obj.a0^(1.0 - 3.0/config.fr_dim);
            obj.amfrac = amfrac_temp * sqrt(0.6);
            obj.bmfrac = 1.0 / config.fr_dim;
            
            % --- Settling constant
            obj.setcon = (2.0/9.0) * obj.del_rho / config.rho_fl * ...
                         config.g / config.kvisc;

            % ============================================================
            % === Step 2 additions: viscosity (T,S) and gamma (ε,ν) ===
            % ============================================================

            % Optional update of viscosity based on temperature & salinity
            if isprop(config, 'salt') && isprop(config, 'temp')
                % Convert temperature to °C if needed
                if config.temp > 200
                    T = config.temp - 273.15;
                else
                    T = config.temp;
                end
                S = config.salt;

                % Dynamic viscosity of pure water (Pa·s)
                mu_w = (2.414e-5) * 10^(247.8 / (T + 133.15));
                % Simple salinity correction
                mu = mu_w * (1 + 0.002 * (S - 35));

                % Convert to kinematic viscosity [cm^2 s^-1]
                rho = 1025;                   % seawater density [kg/m^3]
                nu_m2s = mu / rho;            % [m^2/s]
                config.kvisc = nu_m2s * 1e4;  % [cm^2/s]
            end

            % Compute shear rate gamma from epsilon and viscosity
            if isprop(config, 'epsilon')
                nu_m2s = config.kvisc * 1e-4;      % cm^2/s → m^2/s
                config.gamma = sqrt(config.epsilon / nu_m2s); % [s^-1]
            end

            % Print for verification
            fprintf('DerivedGrid → viscosity = %.3e cm^2/s, gamma = %.3e s^-1\n', ...
                config.kvisc, config.gamma);
        end
        
        % ============================================================
        % === Existing radius/diameter helper methods (unchanged) ===
        % ============================================================

        function r_i = getFractalRadii(obj)
            r_i = obj.amfrac * obj.av_vol.^obj.bmfrac;
        end
        
        function r_v = getConservedRadii(obj)
            r_v = (0.75/pi * obj.av_vol).^(1.0/3.0);
        end
        
        function d_i = getImageDiameters(obj, config)
            r_i = obj.getFractalRadii();
            d_i = 2.0 * config.r_to_rg * r_i;
        end
        
        function d_v = getVolumeDiameters(obj)
            r_v = obj.getConservedRadii();
            d_v = 2.0 * r_v;
        end
    end
end