classdef SettlingVelocityService < handle
    %SETTLINGVELOCITYSERVICE Service for calculating particle settling velocities
    
    methods (Static)
        function v = velocity(r, rcons, setcon)
            %VELOCITY Calculate settling velocities of particles
            % v = particle settling velocities [cm s^{-1}]
            % r = particle radii [cm]
            % rcons = radii of particle conserved volumes [cm]
            % setcon = (2g/9eta)*(delta_rho/rho_fluid)
            
            v = setcon * rcons .* rcons .* rcons ./ r;
        end
        
        function v = velocityWithConfig(r, rcons, config, grid)
            %VELOCITYWITHCONFIG Calculate settling velocities using config and grid
            % Convenience method that extracts setcon from grid
            
            v = SettlingVelocityService.velocity(r, rcons, grid.setcon);
        end
        
        function v = velocityForSections(grid)
            %VELOCITYFORSECTIONS Calculate settling velocities for all sections
            % Returns settling velocities for the entire grid
            
            r_i = grid.getFractalRadii();
            r_v = grid.getConservedRadii();
            
            v = SettlingVelocityService.velocity(r_i, r_v, grid.setcon);
        end
    end
end
