classdef LinearProcessBuilder < handle
    %LINEARPROCESSBUILDER Builds linear operators for growth, sinking, and disaggregation

    methods (Static)
        function G = growthMatrix(config, grid)
            %GROWTHMATRIX Build growth matrix for section-wise growth transfers
            % Returns banded matrix for growth transfers

            n_sections = config.n_sections;
            growth_loss = zeros(n_sections, 1);
            growth_gain = zeros(n_sections-1, 1);

            if config.gro_sec > 0
                growth_loss(config.gro_sec:n_sections-1) = -1;
                growth_gain(config.gro_sec:end) = 2;
            end

            G = diag(growth_loss) + diag(growth_gain, -1);
            G(1,1) = 1;

            G = config.growth * G;
        end

        function S = sinkingMatrix(config, grid)
            %SINKINGMATRIX Build sinking loss matrix
            % Returns diagonal matrix for settling losses

            fractal_radius = grid.getFractalRadii();
            conserved_radius = grid.getConservedRadii();

            settling_vel = SettlingVelocityService.velocity(fractal_radius, conserved_radius, grid.setcon);
            settling_vel = settling_vel / 100 * config.day_to_sec / config.dz;

            S = diag(settling_vel);
        end

        function [Dminus, Dplus] = disaggregationMatrices(config)
            %DISAGGREGATIONMATRICES Build disaggregation matrices
            % Returns Dminus (loss) and Dplus (gain) matrices

            n_sections = config.n_sections;

            % Match legacy CalcCoagDeriv loop:
            % for isec = 2 : n_sections-1
            %   dvdt(isec) = dvdt(isec) - c3*c4^isec*(v(isec) - c4*v(isec+1));
            % end
            % So only sections 2..n-1 participate
            Dminus = zeros(n_sections);
            Dplus = zeros(n_sections);
            if n_sections > 2
                idx = 2:(n_sections-1);
                for k = idx
                    Dminus(k, k) = config.c3 * config.c4^k;
                    Dplus(k, k-1) = config.c3 * config.c4^(k+1);
                end
            end
        end

        function L = linearMatrix(config, grid)
            %LINEARMATRIX Build combined linear matrix (growth - sinking)
            % Returns G - S matrix

            G = LinearProcessBuilder.growthMatrix(config, grid);
            S = LinearProcessBuilder.sinkingMatrix(config, grid);

            L = G - S;
        end
    end
end
