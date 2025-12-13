classdef LinearProcessBuilder < handle
    %LINEARPROCESSBUILDER Builds operators for Growth, Sinking (1D), and Breakup
    
    methods (Static)
        function G = growthMatrix(config, grid)
            % Standard 0-D Growth (Internal to a section, assumes constant across depth later)
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
            % SINKINGMATRIX (1-D Transport)
            % Builds a Global Sparse Matrix for Vertical Sinking.
            % Moves mass from Layer k -> Layer k+1
            
            % 1. Calculate Sinking Velocity w_s [cm/s] per section
            fractal_radius = grid.getFractalRadii();
            conserved_radius = grid.getConservedRadii();
            w_s_cms = SettlingVelocityService.velocity(fractal_radius, conserved_radius, grid.setcon);
            
            % 2. Convert to Transport Rate [1/day]
            % Rate = Velocity / Distance
            % CAREFUL: dz is in meters, w_s is in cm/s
            dz_cm = config.dz * 100;
            rate_per_sec = w_s_cms / dz_cm;
            rate_per_day = rate_per_sec * config.day_to_sec; % [d^-1]
            
            % 3. Define Grid Dimensions
            Nz = floor(config.z_max / config.dz);
            Ns = config.n_sections;
            N_total = Nz * Ns;
            
            if ~config.use_column || Nz <= 1
                % Fallback for 0-D Slab: Just loss term (Sink out of box)
                % Returns simple Ns x Ns diagonal
                S = diag(rate_per_day);
                return;
            end
            
            fprintf('  > Building 1D Transport Matrix (%d layers, %d bins)...\n', Nz, Ns);

            % 4. Build Sparse Matrix Inter-Layer Connections
            % We order the state vector as: [Layer1(all bins); Layer2(all bins); ...]
            
            % --- A. Loss from current layer (Diagonal) ---
            % All layers (1 to Nz) lose mass downwards
            diag_loss = repmat(-rate_per_day(:), Nz, 1); 
            S_loss = spdiags(diag_loss, 0, N_total, N_total);
            
            % --- B. Gain from layer above (Lower Block Diagonal) ---
            % Layers 2..Nz gain mass from 1..Nz-1
            % We construct this explicitly to be safe with sparse indices
            
            row_inds = [];
            col_inds = [];
            vals     = [];
            
            % Pre-allocate approximate size for speed (Ns * (Nz-1))
            num_entries = Ns * (Nz - 1);
            row_inds = zeros(num_entries, 1);
            col_inds = zeros(num_entries, 1);
            vals     = zeros(num_entries, 1);
            
            count = 0;
            for z = 2:Nz
                source_z = z - 1;
                % Offset indices
                offset_dest   = (z-1)*Ns;
                offset_source = (source_z-1)*Ns;
                
                for s = 1:Ns
                    count = count + 1;
                    row_inds(count) = offset_dest + s;
                    col_inds(count) = offset_source + s;
                    vals(count)     = rate_per_day(s);
                end
            end
            
            S_gain = sparse(row_inds, col_inds, vals, N_total, N_total);
            
            % Total Transport Operator: Loss + Gain
            S = S_loss + S_gain;
        end

        function [Dminus, Dplus] = disaggregationMatrices(config)
            % Standard 0-D Linear Breakup (Legacy support)
            n_sections = config.n_sections;
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
            % Combines Growth and Sinking into one operator
            G = LinearProcessBuilder.growthMatrix(config, grid);
            S = LinearProcessBuilder.sinkingMatrix(config, grid);
            
            if config.use_column && floor(config.z_max / config.dz) > 1
                % In 1D:
                % S is already Big (Total x Total)
                % G is Small (Ns x Ns) -> Must expand to Big
                
                Nz = floor(config.z_max / config.dz);
                
                % Kronecker product: Repeat G for every layer
                G_big = kron(speye(Nz), sparse(G)); 
                
                % Combine: L = Growth(at every z) + Transport(between z)
                % Note: S already contains the negative signs for loss and positive for gain
                L = G_big + S; 
            else
                % 0-D Case (Slab)
                % S is just a diagonal of positive rates (loss)
                % L = Growth - Loss
                L = G - S; 
            end
        end
    end
end