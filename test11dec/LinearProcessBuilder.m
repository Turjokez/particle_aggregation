classdef LinearProcessBuilder < handle
    %LINEARPROCESSBUILDER Builds linear operators for growth, sinking, and disaggregation
    %
    % ---------------------------------------------------------------------
    % 
    % NEW-2025-12-11: Still operates on a 0-D slab state vector of length
    % n_sections. The vertical (1-D column) structure is carried by
    % SimulationConfig / DerivedGrid but is not yet expanded into the
    % operator dimensions here. That means:
    %   - Growth is purely sectional.
    %   - Sinking is a bulk loss from the slab with an effective dz.
    % A full 1-D (Ns * Nz) operator would require a larger re-write.
    % ---------------------------------------------------------------------
    %
    % NEW-2025-12-11 (UPDATED):
    % This file now supports BOTH:
    %   - 0-D slab: state size = Ns
    %   - 1-D column: state size = Ns*Nz, stored layer-by-layer as:
    %       v = [N(:,1); N(:,2); ...; N(:,Nz)]  (MATLAB column-major on Ns x Nz)
    %
    % For 1-D column:
    %   - growth acts only in the surface layer (layer 1)
    %   - sinking is vertical transfer between layers + export out bottom
    %
    % NEW-2025-12-12:
    % Add optional config.growth_mode:
    %   'shift' (legacy) = your old section-transfer growth using gro_sec
    %   'pp'             = primary production source: +mu*N1 only
    % If growth_mode is not present, defaults to 'shift' to preserve legacy.
    %
    % NEW-2025-12-12 (UPDATED AGAIN):
    % Add optional config.sinking_form:
    %   'loop' (legacy)  = your for-loop block assembly
    %   'flux' (new)     = conservative flux-form operator (recommended)
    % Default = 'flux' for clarity, but we keep loop version fully.
    %
    % NEW-2025-12-14 (Adrian-style PP):
    % In growth_mode='pp', PP is added explicitly in CoagulationRHS (mu*N1_surface).
    % So growthMatrix should return ZERO in pp mode (avoid double counting).

    methods (Static)

        function G = growthMatrix(config, grid)
            %GROWTHMATRIX Build growth matrix for section-wise growth transfers
            % Returns banded matrix for growth transfers
            %
            % NEW-2025-12-11:
            %   If use_column = true, returns (Ns*Nz x Ns*Nz) where growth
            %   acts only in the first (surface) layer and is zero below.
            %
            % NEW-2025-12-12:
            %   Optional config.growth_mode = 'pp' to support Adrian tests:
            %     Input = mu * N1(t)
            %   Implemented as matrix with only G(1,1)=mu (surface only if column).

            n_sections = config.n_sections;

            % ==========================================================
            % NEW-2025-12-12: decide growth mode (default = legacy 'shift')
            % ==========================================================
            growth_mode = 'shift';
            if isprop(config, 'growth_mode') && ~isempty(config.growth_mode)
                growth_mode = lower(string(config.growth_mode));
            end

            % -------------------------
            % Original 0-D growth block
            % -------------------------
            growth_loss = zeros(n_sections, 1);
            growth_gain = zeros(n_sections-1, 1);

            if config.gro_sec > 0
                growth_loss(config.gro_sec:n_sections-1) = -1;
                growth_gain(config.gro_sec:end) = 2;
            end

            G0_shift = diag(growth_loss) + diag(growth_gain, -1);

            % ----------------------------------------------------------
            % OLD (kept):
            % G0(1,1) = 1;
            % G0 = config.growth * G0;
            % ----------------------------------------------------------

            % NEW (explicit naming):
            G0_shift(1,1) = 1;                 % legacy boundary trick (kept)
            G0_shift = config.growth * G0_shift;

            % ==========================================================
            % NEW-2025-12-12: PP-style growth matrix (Adrian Test 1)
            % ==========================================================
            % Primary production acts like exponential source on first size class:
            %   dN1/dt += mu * N1
            % No transfers to other sections from "growth" operator here.
            G0_pp = sparse(n_sections, n_sections);

            % -----------------------------
            % OLD (kept but disabled):
            % -----------------------------
            % G0_pp(1,1) = config.growth;   % mu on first section only
            %
            % -----------------------------
            % NEW-2025-12-14:
            % PP is handled in CoagulationRHS (mu*N1_surface), so keep this ZERO.
            % -----------------------------
            G0_pp(1,1) = 0;

            % Pick which base matrix to use in 0-D
            if growth_mode == "pp"
                G0 = G0_pp;
            else
                % OLD behavior by default
                % OLD (kept conceptually):
                % G0 = config.growth * (diag(growth_loss) + diag(growth_gain, -1));
                G0 = G0_shift;
            end

            % ----------------------------------------------------------
            % NEW-2025-12-11: Column version (surface-only growth)
            % ----------------------------------------------------------
            if ~isprop(config, 'use_column') || ~config.use_column
                G = G0;  % original behavior (or PP mode)
                return;
            end

            Nz = config.getNumLayers();

            % Build block matrix: first layer gets G0, deeper layers get 0
            % v ordering is [layer1; layer2; ...] with each layer size Ns
            G = sparse(n_sections*Nz, n_sections*Nz);
            G(1:n_sections, 1:n_sections) = sparse(G0);
        end


        function [S, sink_rate, settling_vel_cmday] = sinkingMatrix(config, grid)
            %SINKINGMATRIX Build sinking loss / transfer matrix
            %
            % Original slab behaviour:
            %   settling_vel [cm s^-1] -> [cm d^-1] -> divide by dz [cm]
            %   S_j = w_s / dz, acting as a first-order loss from the slab.
            %
            % NEW-2025-12-11:
            %   If use_column = true, build a vertical transfer operator:
            %     dN_k/dt includes:
            %       - (w/dz) * N_k          (loss from layer k)
            %       + (w/dz) * N_{k-1}      (gain from layer above)
            %   Bottom layer loss leaves the system (export).
            %
            % NEW-2025-12-11 (small addition):
            %   Also return sink_rate (Ns x 1, in 1/day), so other files can
            %   compute bottom export flux cleanly without reading matrix diagonals.
            %
            % NEW-2025-12-12:
            %   Added optional third output: settling_vel_cmday (Ns x 1).
            %   This helps you debug whether you should use:
            %     export flux  F = N_bottom .* w   (cm/day form)
            %   vs
            %     export loss  L = (w/dz) .* N_bottom  (1/day form)

            fractal_radius   = grid.getFractalRadii();
            conserved_radius = grid.getConservedRadii();

            settling_vel = SettlingVelocityService.velocity(fractal_radius, conserved_radius, grid.setcon);
            % cm/s -> cm/day
            settling_vel_cmday = settling_vel * config.day_to_sec;

            % ==========================================================
            % NEW-2025-12-11: interpret dz as layer thickness [m]
            % convert dz from m to cm
            % ==========================================================
            dz_eff = config.dz;          % [m]
            dz_eff_cm = dz_eff * 100;    % [cm]

            % Convert to first-order rate [1/day]
            sink_rate = settling_vel_cmday / dz_eff_cm;   % [1/day], size Ns x 1

            % -------------------------
            % 0-D slab (original)
            % -------------------------
            if ~isprop(config, 'use_column') || ~config.use_column
                S = diag(sink_rate);    % positive diagonal (loss)
                return;
            end

            % ----------------------------------------------------------
            % NEW-2025-12-12: choose sinking form
            % ----------------------------------------------------------
            sinking_form = "flux";
            if isprop(config, 'sinking_form') && ~isempty(config.sinking_form)
                sinking_form = lower(string(config.sinking_form));
            end

            Nz = config.getNumLayers();
            Ns = config.n_sections;

            W = spdiags(sink_rate, 0, Ns, Ns);  % positive rates on diagonal
            Ntot = Ns * Nz;

            % ----------------------------------------------------------
            % OLD (kept): loop-built block matrix
            % ----------------------------------------------------------
            % IMPORTANT CLARIFICATION (NEW-2025-12-11):
            % We build S so that L = G - S produces:
            %   dN1/dt += -W*N1
            %   dNk/dt += +W*N_{k-1} - W*N_k   (k>1)
            %
            % That means:
            %   - diagonal blocks of S are +W
            %   - subdiagonal blocks of S are -W
            %
            % OLD LOOP FORM:
            % S = sparse(Ntot, Ntot);
            % for k = 1:Nz
            %     rows = (k-1)*Ns + (1:Ns);
            %     cols = rows;
            %
            %     % +W on diagonal block
            %     S(rows, cols) = S(rows, cols) + W;
            %
            %     % -W on subdiagonal block
            %     if k >= 2
            %         cols_above = (k-2)*Ns + (1:Ns);
            %         S(rows, cols_above) = S(rows, cols_above) - W;
            %     end
            % end

            if sinking_form == "loop"
                % NEW-2025-12-12: keep your loop as an option (same math)
                S = sparse(Ntot, Ntot);

                for k = 1:Nz
                    rows = (k-1)*Ns + (1:Ns);
                    cols = rows;

                    % +W on diagonal block
                    S(rows, cols) = S(rows, cols) + W;

                    % -W on subdiagonal block
                    if k >= 2
                        cols_above = (k-2)*Ns + (1:Ns);
                        S(rows, cols_above) = S(rows, cols_above) - W;
                    end
                end
            else
                % ------------------------------------------------------
                % NEW-2025-12-12: conservative flux-form operator (recommended)
                % ------------------------------------------------------
                Iz = speye(Nz);
                shift_down = spdiags(ones(Nz,1), -1, Nz, Nz);  % ones on subdiagonal
                Tz = Iz - shift_down;

                S = kron(Tz, W);
            end

            % ----------------------------------------------------------
            % Optional debug check (does not change results)
            % Turn on by setting config.debug_sinking = true
            % ----------------------------------------------------------
            if isprop(config, 'debug_sinking') && config.debug_sinking
                vtest = abs(rand(Ntot,1));
                dtest = -(S * vtest);  % because L includes -S part

                Nmat = reshape(vtest, [Ns, Nz]);
                Nbottom = Nmat(:, end);
                export_loss_rate = sum(sink_rate .* Nbottom);

                colsum = sum(dtest);
                fprintf('[debug_sinking] sum(dN/dt) = %.6e,  -export_loss_rate = %.6e\n', colsum, -export_loss_rate);
            end
        end


        function [Dminus, Dplus] = disaggregationMatrices(config)
            %DISAGGREGATIONMATRICES Build disaggregation matrices
            % Returns Dminus (loss) and Dplus (gain) matrices
            %
            % IMPORTANT NEW-2025-12-12:
            % Your legacy matrix uses a SUBDIAGONAL gain (k,k-1).
            % But your OLD manual disagg loop (commented in CoagulationRHS) was:
            %   dv(k) += -a*v(k) + a*c4*v(k+1)
            % which corresponds to a SUPERDIAGONAL gain (k,k+1).
            %
            % To avoid breaking old behavior, we DO NOT delete anything.
            % We add a switch:
            %   config.disagg_mode = 'legacy'    (default if missing)
            %                      = 'loop_equiv'  (matches the old loop form)

            n_sections = config.n_sections;

            % ----------------------------------------------------------
            % NEW: choose disaggregation mode
            % ----------------------------------------------------------
            disagg_mode = "legacy";  % default keeps your current behavior
            if isprop(config,'disagg_mode') && ~isempty(config.disagg_mode)
                disagg_mode = lower(string(config.disagg_mode));
            end

            % ----------------------------------------------------------
            % Allocate
            % ----------------------------------------------------------
            Dminus = zeros(n_sections);
            Dplus  = zeros(n_sections);

            % ==========================================================
            % MODE 1: LEGACY (YOUR CODE — UNCHANGED)
            % ==========================================================
            if disagg_mode == "legacy"
                % --------- BEGIN ORIGINAL LEGACY CODE (DO NOT CHANGE) ----------
                Dminus = zeros(n_sections);
                Dplus = zeros(n_sections);
                if n_sections > 2
                    idx = 2:(n_sections-1);
                    for k = idx
                        Dminus(k, k)   = config.c3 * config.c4^k;
                        Dplus(k, k-1)  = config.c3 * config.c4^(k+1);
                    end
                end
                % --------- END ORIGINAL LEGACY CODE (DO NOT CHANGE) ------------

                % NEW-2025-12-11:
                % NOTE: For column runs, these remain sectional (Ns x Ns) and
                % CoagulationRHS will apply them layer-by-layer.
                return;
            end

            % ==========================================================
            % MODE 2: LOOP-EQUIVALENT (matches old manual loop form)
            % dv(k) += -a*v(k) + a*c4*v(k+1)
            % where a = c3*c4^k, for k = 2..Ns-1
            % ==========================================================
            % This fixes the "wrong diagonal direction" issue and prevents
            % disaggregation from acting like an unintended source.
            %
            % NOTE:
            % We intentionally keep the same a(k)=c3*c4^k structure you used.
            % Only the destination index is corrected to be consistent with:
            %   (vk(isec) - c4*vk(isec+1)) term
            if n_sections > 2
                idx = 2:(n_sections-1);
                for k = idx
                    a = config.c3 * config.c4^k;
                    Dminus(k, k)   = a;              % -a*v(k)
                    Dplus(k, k+1)  = a * config.c4;  % +a*c4*v(k+1)
                end
            end

            % Optional: use sparse (safe, but not required)
            % Dminus = sparse(Dminus);
            % Dplus  = sparse(Dplus);

            % NEW-2025-12-11:
            % NOTE: For column runs, these remain sectional (Ns x Ns) and
            % CoagulationRHS will apply them layer-by-layer.
        end


        function L = linearMatrix(config, grid)
            %LINEARMATRIX Build combined linear matrix (growth - sinking)
            % Returns G - S matrix

            G = LinearProcessBuilder.growthMatrix(config, grid);

            % OLD:
            % S = LinearProcessBuilder.sinkingMatrix(config, grid);

            % NEW-2025-12-11: also capture sink_rate if needed later
            % [S, ~] = LinearProcessBuilder.sinkingMatrix(config, grid);  % NEW-2025-12-11

            % NEW-2025-12-12: keep extra output for debugging if needed
            [S, ~, ~] = LinearProcessBuilder.sinkingMatrix(config, grid);  % NEW-2025-12-12

            L = G - S;
        end
    end
end