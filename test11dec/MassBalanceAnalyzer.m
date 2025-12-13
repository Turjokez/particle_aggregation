classdef MassBalanceAnalyzer < handle
    %MASSBALANCEANALYZER Analyzes mass balance for coagulation simulation
    %
    % NEW-2025-12-11:
    % Supports column states where spec is [time x (Ns*Nz)].
    % In column mode:
    %   - internal sinking between layers is NOT treated as loss
    %   - only export out the bottom layer is counted as settling loss

    methods (Static)
        function [gains, losses] = sectional(spec, operators)
            %SECTIONAL Calculate sectional mass balance
            % spec = concentration matrix (time x sections) OR (time x Ns*Nz)
            % operators = struct with coagulation and linear operators
            % Returns: gains, losses structs with sectional data
            % NOTE: This replicates the exact legacy algorithm including bugs
            %
            % NEW-2025-12-11:
            % Column detection by (operators.config) or by dimension fallback.

            [n_times, n_cols] = size(spec);

            % ----------------------------------------------------------
            % NEW-2025-12-11: Column detection
            % ----------------------------------------------------------
            use_column = false;
            Ns = [];
            Nz = [];

            if isfield(operators, 'config') && isprop(operators.config, 'use_column') && operators.config.use_column
                use_column = true;
                Ns = operators.config.n_sections;
                Nz = operators.config.getNumLayers();
            else
                % Fallback detection by dimension if config not provided
                if isfield(operators, 'betas')
                    Ns_guess = operators.betas.getNumSections();
                    if mod(n_cols, Ns_guess) == 0 && n_cols ~= Ns_guess
                        use_column = true;
                        Ns = Ns_guess;
                        Nz = n_cols / Ns_guess;
                    end
                end
            end

            if use_column
                % ======================================================
                % NEW-2025-12-11: Column branch
                % ======================================================

                coag_gains  = zeros(n_times, Ns);
                coag_losses = zeros(n_times, Ns);

                growth_gain = zeros(n_times, Ns);
                growth_loss = zeros(n_times, Ns);

                sink_losses = zeros(n_times, Ns);

                % ------------------------------------------------------
                % NEW-2025-12-11: Extract Ns×Ns surface blocks when the
                % operators are expanded to Ns*Nz × Ns*Nz
                % ------------------------------------------------------
                Gop = operators.growth;
                Sop = operators.sink_loss;

                % OLD (kept): assumes Ns×Ns
                % G = operators.growth;
                % sink_rate = diag(operators.sink_loss);

                % NEW-2025-12-11: robust block extraction
                if size(Gop,1) == Ns*Nz
                    G0 = Gop(1:Ns, 1:Ns);        % surface layer growth block
                else
                    G0 = Gop;                    % already Ns×Ns (older runs)
                end

                if size(Sop,1) == Ns*Nz
                    W0 = Sop(1:Ns, 1:Ns);        % this is +W (rates) on surface block
                    sink_rate = diag(W0);        % Ns×1
                else
                    sink_rate = diag(Sop);       % Ns×1 in slab mode
                end
                
                % NOTE:
                % In column mode, every diagonal block of Sop contains the same +W (rates = w/dz),
                % so reading from the surface block is OK and equals the bottom block numerically.

                % Growth diag vectors (legacy-style)
                g_gain_vec = diag(G0, -1);
                g_gain_vec = [G0(1,1); g_gain_vec];  % size Ns
                g_loss_vec = diag(G0);
                g_loss_vec(1) = 0.0;

                for it = 1:n_times
                    vflat = spec(it, :)';

                    % reshape to Ns x Nz (same ordering as model)
                    N2 = reshape(vflat, [Ns, Nz]);     % Ns x Nz

                    % -------- coagulation per layer, then sum over z --------
                    cg = zeros(Ns, 1);
                    cl = zeros(Ns, 1);

                    for k = 1:Nz
                        vcon_r = N2(:, k)';              % 1 x Ns
                        vcon_shift = [0, vcon_r(1:Ns-1)];

                        % Gains (legacy-like)
                        term1 = vcon_r * operators.betas.b2;
                        term1 = vcon_r .* term1;

                        term2 = vcon_r * operators.betas.b1;
                        term2 = term2 .* vcon_shift;

                        cg = cg + (term1 + term2)';

                        % Losses (legacy-like)
                        term3 = vcon_r * (operators.betas.b3 + operators.betas.b4 + operators.betas.b5);
                        term3 = vcon_r .* term3;

                        cl = cl + term3';
                    end

                    coag_gains(it, :)  = cg';
                    coag_losses(it, :) = cl';

                    % -------- growth diagnostics on column-integrated spec --------
                    vcol = sum(N2, 2);  % Ns x 1

                    % OLD NOTE (kept): slab version replicates a bug
                    % In column mode we do NOT replicate that bug.
                    growth_gain(it, :) = (g_gain_vec .* vcol)';
                    growth_loss(it, :) = (g_loss_vec .* vcol)';

                % -------- sinking export loss: only bottom layer --------
                vbot = N2(:, end);  % Ns x 1
                
                % OLD (kept): this is a RATE-like loss (1/day * conc), not per-area flux
                % sink_losses(it, :) = (sink_rate .* vbot)';
                
                % NEW-2025-12-12: convert to bottom export FLUX per area [cm^3 m^-2 d^-1]
                % sink_rate = w/dz [1/day], so F = (w/dz)*N * dz * 1e4  -> dz cancels if sink_rate already includes /dz
                sink_losses(it, :) = (1e4 * (sink_rate .* vbot))';

                end

                gains.coag = coag_gains;
                gains.growth = growth_gain;

                losses.coag = coag_losses;
                losses.growth = growth_loss;
                losses.settl = sink_losses;

                % Optional metadata (safe)
                gains.Ns = Ns; %#ok<STRNU>
                gains.Nz = Nz; %#ok<STRNU>
                return;
            end

            % ==========================================================
            % ORIGINAL 0-D slab code (unchanged below)
            % ==========================================================
            [n_times, n_sections] = size(spec);

            coag_gains = zeros(n_times, n_sections);
            coag_losses = zeros(n_times, n_sections);

            for i_time = 1:n_times
                vcon_r = spec(i_time, :);
                vcon_shift = [0, vcon_r(1:n_sections-1)];

                term1 = vcon_r * operators.betas.b2;
                term1 = vcon_r .* term1;

                term2 = vcon_r * operators.betas.b1;
                term2 = term2 .* vcon_shift;

                coag_gains(i_time, :) = term1 + term2;

                term3 = vcon_r * (operators.betas.b3 + operators.betas.b4 + operators.betas.b5);
                term3 = vcon_r .* term3;

                coag_losses(i_time, :) = term3;
            end

            sinking = diag(operators.sink_loss);
            sinking = sinking';
            sinking = sinking(ones(n_times, 1), :);
            sink_losses = sinking .* spec;

            g_gain = diag(operators.growth, -1);
            g_gain = [operators.growth(1, 1); g_gain];
            g_gain = g_gain';
            g_gain = g_gain(ones(i_time, 1), :);  % BUG: legacy behavior

            g_loss = diag(operators.growth);
            g_loss(1) = 0.0;
            g_loss = g_loss';
            g_loss = g_loss(ones(i_time, 1), :);  % BUG: legacy behavior

            growth_gain = g_gain .* spec;
            growth_loss = g_loss .* spec;

            gains.coag = coag_gains;
            gains.growth = growth_gain;

            losses.coag = coag_losses;
            losses.growth = growth_loss;
            losses.settl = sink_losses;
        end

        function [total_gains, total_losses] = total(spec, operators)
            %TOTAL Calculate total mass balance
            %
            % NEW-2025-12-11:
            % In column mode, settling losses = export out bottom only.

            [n_times, n_cols] = size(spec);

            use_column = false;
            Ns = [];
            Nz = [];

            if isfield(operators, 'config') && isprop(operators.config, 'use_column') && operators.config.use_column
                use_column = true;
                Ns = operators.config.n_sections;
                Nz = operators.config.getNumLayers();
            else
                if isfield(operators, 'betas')
                    Ns_guess = operators.betas.getNumSections();
                    if mod(n_cols, Ns_guess) == 0 && n_cols ~= Ns_guess
                        use_column = true;
                        Ns = Ns_guess;
                        Nz = n_cols / Ns_guess;
                    end
                end
            end

            if use_column
                % ======================================================
                % NEW-2025-12-11: Column totals
                % ======================================================

                % Extract sink rates from surface block if expanded
                Sop = operators.sink_loss;
                if size(Sop,1) == Ns*Nz
                    W0 = Sop(1:Ns, 1:Ns);
                    sink_rate = diag(W0);
                else
                    sink_rate = diag(Sop);
                end

                % Extract growth surface block if expanded
                Gop = operators.growth;
                if size(Gop,1) == Ns*Nz
                    G0 = Gop(1:Ns, 1:Ns);
                else
                    G0 = Gop;
                end

                total_sink_losses = zeros(n_times, 1);
                net_growth = zeros(n_times, 1);
                coag_losses = zeros(n_times, 1);

                for it = 1:n_times
                    vflat = spec(it, :)';
                    N2 = reshape(vflat, [Ns, Nz]);   % Ns x Nz

                    % export out bottom only
                    vbot = N2(:, end);
                    
                    % OLD (kept): rate-form (not per-area flux)
                    % total_sink_losses(it) = sum(sink_rate .* vbot);
                    
                    % NEW-2025-12-12: bottom export flux per area [cm^3 m^-2 d^-1]
                    total_sink_losses(it) = 1e4 * sum(sink_rate .* vbot);

                    % net growth: apply growth operator to column-integrated spec
                    vcol = sum(N2, 2);
                    g1 = G0 * vcol;
                    net_growth(it) = sum(g1);

                    % coagulation "loss from system" diagnostic (legacy-style)
                    v = vcol(:)';
                    n_sections = Ns;

                    term1 = operators.betas.b4(n_sections, n_sections) * v(n_sections) * v(n_sections);

                    term2 = (operators.betas.b5(n_sections, :) .* v) * v(n_sections);
                    term2 = term2';

                    term3 = (operators.betas.b2(:, n_sections) .* v') * v(n_sections);

                    term4 = term2 - term3;
                    term4 = sum(term4');

                    term5 = (operators.betas.b3(:, n_sections) .* v') * v(n_sections);
                    term5 = sum(term5');

                    coag_losses(it) = term1 + term4 + term5;
                end

                total_gains.growth = net_growth;
                total_losses.sett = total_sink_losses;
                total_losses.coag = coag_losses;
                return;
            end

            % ==========================================================
            % ORIGINAL 0-D slab code (unchanged below)
            % ==========================================================
            [n_times, n_sections] = size(spec);

            sinking = diag(operators.sink_loss);
            sinking = sinking';
            sinking = sinking(ones(n_times, 1), :);
            sink_losses = sinking .* spec;
            total_sink_losses = sum(sink_losses, 2);

            net_growth = zeros(n_times, 1);
            for i_time = 1:n_times
                v = spec(i_time, :)';
                g1 = operators.growth * v;
                net_growth(i_time) = sum(g1);
            end

            coag_losses = zeros(n_times, 1);
            for i_time = 1:n_times
                v = spec(i_time, :);
                v_r = v';

                term1 = operators.betas.b4(n_sections, n_sections) * v(n_sections) * v(n_sections);

                term2 = (operators.betas.b5(n_sections, :) .* v) * v(n_sections);
                term2 = term2';

                term3 = (operators.betas.b2(:, n_sections) .* v') * v(n_sections);

                term4 = term2 - term3;
                term4 = sum(term4');

                term5 = (operators.betas.b3(:, n_sections) .* v') * v(n_sections);
                term5 = sum(term5');

                coag_losses(i_time) = term1 + term4 + term5;
            end

            total_gains.growth = net_growth;
            total_losses.sett = total_sink_losses;
            total_losses.coag = coag_losses;
        end

        function [gains, losses] = sectionalWithRates(spec, rhs)
            %SECTIONALWITHRATES Calculate sectional mass balance using RHS

            [n_times, n_sections] = size(spec);

            coag_gains = zeros(n_times, n_sections);
            coag_losses = zeros(n_times, n_sections);

            for i_time = 1:n_times
                [term1, term2, term3, term4, term5] = rhs.rateTerms(spec(i_time, :)');

                coag_gains(i_time, :) = term1 + term2;
                coag_losses(i_time, :) = -term4;
            end

            gains.coag = coag_gains;
            gains.growth = zeros(n_times, n_sections);

            losses.coag = coag_losses;
            losses.growth = zeros(n_times, n_sections);
            losses.settl = zeros(n_times, n_sections);
        end

        function displayBalanceSummary(gains, losses, time_vector)
            %DISPLAYBALANCESUMMARY Display summary of mass balance

            fprintf('Mass Balance Summary:\n');
            fprintf('  Time points: %d\n', length(time_vector));
            fprintf('  Time range: [%.2f, %.2f]\n', time_vector(1), time_vector(end));

            if isfield(gains, 'coag')
                fprintf('  Coagulation gains range: [%.2e, %.2e]\n', ...
                    min(gains.coag(:)), max(gains.coag(:)));
            end
            if isfield(gains, 'growth')
                fprintf('  Growth gains range: [%.2e, %.2e]\n', ...
                    min(gains.growth(:)), max(gains.growth(:)));
            end

            if isfield(losses, 'coag')
                fprintf('  Coagulation losses range: [%.2e, %.2e]\n', ...
                    min(losses.coag(:)), max(losses.coag(:)));
            end
            if isfield(losses, 'settl')
                fprintf('  Settling losses range: [%.2e, %.2e]\n', ...
                    min(losses.settl(:)), max(losses.settl(:)));
            end
        end
    end
end
