classdef MassBalanceAnalyzer < handle
    %MASSBALANCEANALYZER Analyzes mass balance for coagulation simulation
    %
    % NEW-2025-12-11:
    % Supports column states where spec is [time x (Ns*Nz)].
    % In column mode:
    %   - internal sinking between layers is NOT treated as loss
    %   - only export out the bottom layer is counted as settling loss
    %
    % NEW-2025-12-21 (FIX):
    % If sink_loss operator diagonal is a RATE w/dz [1/day],
    % then export flux per area is:
    %   F = (w/dz)*N_bottom * dz * 1e4  = w*N_bottom*1e4
    % so we MUST multiply by dz_cm*1e4 when using sink_rate from operator.
    %
    % NEW-2025-12-21 (UPDATED FIX):
    % Prefer computing bottom export from settling velocity (cm/day):
    %   F_bottom = 1e4 * sum( N_bottom .* w )
    % This matches your verified check (F_recalc vs F_stored).
    %
    % NEW-2026-01-20:
    % Add consistent bin-weighting for export diagnostics in column mode:
    %   export_weight="ones" or "vbin"
    % and compute export as:
    %   Fbin = 1e4*(w_cmday .* vbot .* wbin)

    methods (Static)

        % ==============================================================
        % NEW-2025-12-21: helper to get settling velocity w (cm/day)
        % ==============================================================
        function w = getSettlingVelCmDay(operators, Ns)
            % Returns Ns x 1 w in cm/day if present, else [].
            w = [];
            if isfield(operators,'settling_vel_cmday') && ~isempty(operators.settling_vel_cmday)
                w = operators.settling_vel_cmday(:);
                if numel(w) ~= Ns
                    % keep safe: if size mismatch, ignore
                    w = [];
                end
            elseif isfield(operators,'settling_vel_cmday') && isempty(operators.settling_vel_cmday)
                w = [];
            end
        end

        % ==============================================================
        % NEW-2026-01-20: helper to pick bin weight
        % ==============================================================
        function wbin = getBinWeightFromOperators(operators, Ns)
            % Default: ones
            wbin = ones(Ns,1);

            % Try to read mode from config
            weight_mode = "ones";
            try
                if isfield(operators,'config') && isprop(operators.config,'export_weight') && ~isempty(operators.config.export_weight)
                    weight_mode = string(operators.config.export_weight);
                end
            catch
            end

            if strcmpi(weight_mode,"vbin")
                % Try to get bin volumes from config/grid if available.
                % We keep it safe: if we can't get it, fall back to ones.
                try
                    if isfield(operators,'grid') && ~isempty(operators.grid)
                        gridS = operators.grid;
                        wtmp = [];
                        try
                            if ismethod(gridS,'getBinVolumes')
                                wtmp = gridS.getBinVolumes();
                            end
                        catch
                        end
                        if isempty(wtmp)
                            try
                                wtmp = Disaggregation.getBinVolumes_cm3(gridS);
                            catch
                                wtmp = [];
                            end
                        end
                        if ~isempty(wtmp)
                            wbin = wtmp(:);
                        end
                    end
                catch
                end
            end
        end

        function [gains, losses] = sectional(spec, operators)
            %SECTIONAL Calculate sectional mass balance
            % spec = concentration matrix (time x sections) OR (time x Ns*Nz)
            % operators = struct with coagulation and linear operators
            % Returns: gains, losses structs with sectional data

            [n_times, n_cols] = size(spec);

            % ----------------------------------------------------------
            % Column detection
            % ----------------------------------------------------------
            use_column = false;
            Ns = [];
            Nz = [];
            dz_m  = [];
            dz_cm = [];

            if isfield(operators, 'config') && isprop(operators.config, 'use_column') && operators.config.use_column
                use_column = true;
                Ns = operators.config.n_sections;
                Nz = operators.config.getNumLayers();

                % dz handling
                if isprop(operators.config,'dz')
                    dz_m = operators.config.dz;
                else
                    dz_m = 1; % fallback
                end
                dz_cm = dz_m * 100;
            else
                % Fallback detection by dimension if config not provided
                if isfield(operators, 'betas')
                    Ns_guess = operators.betas.getNumSections();
                    if mod(n_cols, Ns_guess) == 0 && n_cols ~= Ns_guess
                        use_column = true;
                        Ns = Ns_guess;
                        Nz = n_cols / Ns_guess;

                        % dz unknown here -> assume 1m if not provided
                        dz_m  = 1;
                        dz_cm = 100;
                    end
                end
            end

            if use_column
                % ======================================================
                % Column branch
                % ======================================================
                coag_gains  = zeros(n_times, Ns);
                coag_losses = zeros(n_times, Ns);

                growth_gain = zeros(n_times, Ns);
                growth_loss = zeros(n_times, Ns);

                % Settling/export (bottom only) per section
                sink_losses = zeros(n_times, Ns);

                % ------------------------------------------------------
                % Operators (may be expanded to Ns*Nz × Ns*Nz)
                % ------------------------------------------------------
                Gop = operators.growth;
                Sop = operators.sink_loss;

                % Extract surface blocks if expanded
                if size(Gop,1) == Ns*Nz
                    G0 = Gop(1:Ns, 1:Ns);
                else
                    G0 = Gop;
                end

                if size(Sop,1) == Ns*Nz
                    W0 = Sop(1:Ns, 1:Ns);   % typically w/dz [1/day]
                    sink_rate = diag(W0);   % Ns×1
                else
                    sink_rate = diag(Sop);
                end

                % NEW-2025-12-21: get settling velocity (preferred)
                w_cmday = MassBalanceAnalyzer.getSettlingVelCmDay(operators, Ns);
                has_w = ~isempty(w_cmday);

                % NEW-2026-01-20: bin weight
                wbin = ones(Ns,1);
                try
                    % If you can pass grid through operators, do it:
                    % operators.grid = obj.grid in CoagulationSimulation.initializeComponents
                    wbin = MassBalanceAnalyzer.getBinWeightFromOperators(operators, Ns);
                catch
                    wbin = ones(Ns,1);
                end

                % Growth diag vectors (legacy-style)
                g_gain_vec = diag(G0, -1);
                g_gain_vec = [G0(1,1); g_gain_vec];  % Ns×1
                g_loss_vec = diag(G0);
                g_loss_vec(1) = 0.0;

                for it = 1:n_times
                    vflat = spec(it, :)';
                    N2 = reshape(vflat, [Ns, Nz]);   % Ns x Nz

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
                    growth_gain(it, :) = (g_gain_vec .* vcol)';
                    growth_loss(it, :) = (g_loss_vec .* vcol)';

                    % -------- sinking export loss: bottom layer only --------
                    vbot = N2(:, end);  % Ns x 1

                    % NEW-2026-01-20: weighted export (consistent)
                    if has_w
                        sink_losses(it, :) = (1e4 * (w_cmday .* vbot .* wbin))';
                    else
                        % fallback
                        sink_losses(it, :) = (sink_rate .* vbot * dz_cm * 1e4)';
                    end
                end

                gains = struct();
                losses = struct();

                gains.coag   = coag_gains;
                gains.growth = growth_gain;

                losses.coag   = coag_losses;
                losses.growth = growth_loss;
                losses.settl  = sink_losses;   % bottom export (per-section), weighted if w present

                % Optional metadata
                gains.Ns = Ns; %#ok<STRNU>
                gains.Nz = Nz; %#ok<STRNU>
                losses.dz_cm = dz_cm; %#ok<STRNU>
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
            % NEW-2025-12-21 (FIX):
            % In column mode, total settling loss uses dz_cm*1e4 if sink_rate = w/dz.
            %
            % NEW-2025-12-21 (UPDATED FIX):
            % Prefer using settling velocity:
            %   Fbottom = 1e4 * sum(Nbot .* w_cmday)
            %
            % NEW-2026-01-20:
            % If export_weight="vbin", compute:
            %   Fbottom = 1e4 * sum(Nbot .* w_cmday .* wbin)

            [n_times, n_cols] = size(spec);

            use_column = false;
            Ns = [];
            Nz = [];
            dz_m  = [];
            dz_cm = [];

            if isfield(operators, 'config') && isprop(operators.config, 'use_column') && operators.config.use_column
                use_column = true;
                Ns = operators.config.n_sections;
                Nz = operators.config.getNumLayers();

                if isprop(operators.config,'dz')
                    dz_m = operators.config.dz;
                else
                    dz_m = 1;
                end
                dz_cm = dz_m * 100;
            else
                if isfield(operators, 'betas')
                    Ns_guess = operators.betas.getNumSections();
                    if mod(n_cols, Ns_guess) == 0 && n_cols ~= Ns_guess
                        use_column = true;
                        Ns = Ns_guess;
                        Nz = n_cols / Ns_guess;
                        dz_m  = 1;
                        dz_cm = 100;
                    end
                end
            end

            if use_column
                % Extract sink rates (surface block)
                Sop = operators.sink_loss;
                if size(Sop,1) == Ns*Nz
                    W0 = Sop(1:Ns, 1:Ns);
                    sink_rate = diag(W0); % likely w/dz [1/day]
                else
                    sink_rate = diag(Sop);
                end

                % NEW-2025-12-21: settling velocity preferred
                w_cmday = MassBalanceAnalyzer.getSettlingVelCmDay(operators, Ns);
                has_w = ~isempty(w_cmday);

                % NEW-2026-01-20: bin weight
                wbin = ones(Ns,1);
                try
                    wbin = MassBalanceAnalyzer.getBinWeightFromOperators(operators, Ns);
                catch
                    wbin = ones(Ns,1);
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

                    if has_w
                        total_sink_losses(it) = 1e4 * sum(w_cmday .* vbot .* wbin);
                    else
                        total_sink_losses(it) = sum(sink_rate .* vbot) * dz_cm * 1e4;
                    end

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

                total_gains = struct();
                total_losses = struct();

                total_gains.growth = net_growth;
                total_losses.sett  = total_sink_losses;
                total_losses.coag  = coag_losses;
                total_losses.dz_cm = dz_cm; %#ok<STRNU>
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

        function displayBalanceSummary(gains, losses, time_vector)
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

        function closureCheckFromOutput(out)
            %CLOSURECHECKFROMOUTPUT Quick column closure check using output_data
            %
            % Checks dM/dt + Fbottom ~= 0 (should be ~0 ONLY if no sources are on)

            if ~isfield(out,'output_data')
                error('closureCheckFromOutput: out has no output_data');
            end
            if ~isfield(out.output_data,'column_total_inventory_cm3m2')
                error('closureCheckFromOutput: missing column_total_inventory_cm3m2');
            end
            if ~isfield(out.output_data,'bottom_total_flux_cm3m2d')
                error('closureCheckFromOutput: missing bottom_total_flux_cm3m2d');
            end

            M = out.output_data.column_total_inventory_cm3m2(:);
            F = out.output_data.bottom_total_flux_cm3m2d(:);

            % time field name in your out is "time" (not t)
            if isfield(out,'time')
                t = out.time(:);
            else
                t = (0:numel(M)-1)'; % fallback
                warning('closureCheckFromOutput: out.time missing; using index time.');
            end

            dMdt = gradient(M, t);
            res  = dMdt + F;

            figure;
            plot(t, dMdt, t, -F, '--', t, res, '-.');
            xlabel('Time [d]');
            ylabel('cm^3 m^{-2} d^{-1}');
            legend('dM/dt','-F_{bottom}','residual = dM/dt + F');
            title('Column closure check: dM/dt = -F_{bottom} (ONLY if no sources)');

            fprintf('max|res| = %.3e\n', max(abs(res)));
            fprintf('max|res|/max|F| = %.3e\n', max(abs(res))/max(abs(F)));
        end

    end
end