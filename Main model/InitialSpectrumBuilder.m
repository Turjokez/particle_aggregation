classdef InitialSpectrumBuilder < handle
    %INITIALSPECTRUMBUILDER Builds initial particle size spectra
    
    methods (Static)
        function spec_init = initialSpectrum(config, grid)
            %INITIALSPECTRUM Build initial spectrum with equal volume per section
            % Returns initial concentration vector for all sections
            
            spec_init = grid.av_vol(1) * ones(config.n_sections, 1);
            
            tfactor = 10.^(0:config.n_sections-1)';
            spec_init = max(spec_init ./ tfactor, 1.0e-30);
            
            spec_init = spec_init * config.num_1;
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
