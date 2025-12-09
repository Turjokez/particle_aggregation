classdef MassBalanceAnalyzer < handle
    %MASSBALANCEANALYZER Analyzes mass balance for coagulation simulation

    methods (Static)
        function [gains, losses] = sectional(spec, operators)
            %SECTIONAL Calculate sectional mass balance
            % spec = concentration matrix (time x sections)
            % operators = struct with coagulation and linear operators
            % Returns: gains, losses structs with sectional data
            % NOTE: This replicates the exact legacy algorithm including bugs

            [n_times, n_sections] = size(spec);

            % Initialize arrays
            coag_gains = zeros(n_times, n_sections);
            coag_losses = zeros(n_times, n_sections);

            % Calculate coagulation gains and losses for each time
            for i_time = 1:n_times
                vcon_r = spec(i_time, :);
                vcon_shift = [0, vcon_r(1:n_sections-1)];

                % Coagulation gains
                term1 = vcon_r * operators.betas.b2;
                term1 = vcon_r .* term1;

                term2 = vcon_r * operators.betas.b1;
                term2 = term2 .* vcon_shift;

                coag_gains(i_time, :) = term1 + term2;

                % Coagulation losses
                term3 = vcon_r * (operators.betas.b3 + operators.betas.b4 + operators.betas.b5);
                term3 = vcon_r .* term3;

                coag_losses(i_time, :) = term3;
            end

            % Sinking losses
            sinking = diag(operators.sink_loss);
            sinking = sinking';
            sinking = sinking(ones(n_times, 1), :);
            sink_losses = sinking .* spec;

            % Growth gains and losses - EXACT legacy replication including the bug
            % Legacy uses i_time (which equals n_times after the loop) instead of n_times
            g_gain = diag(operators.growth, -1);
            g_gain = [operators.growth(1, 1); g_gain];
            g_gain = g_gain';
            g_gain = g_gain(ones(i_time, 1), :);  % BUG: i_time = n_times here (legacy behavior)

            g_loss = diag(operators.growth);
            g_loss(1) = 0.0;
            g_loss = g_loss';
            g_loss = g_loss(ones(i_time, 1), :);  % BUG: i_time = n_times here (legacy behavior)

            growth_gain = g_gain .* spec;
            growth_loss = g_loss .* spec;

            % Assemble results
            gains.coag = coag_gains;
            gains.growth = growth_gain;

            losses.coag = coag_losses;
            losses.growth = growth_loss;
            losses.settl = sink_losses;
        end

        function [total_gains, total_losses] = total(spec, operators)
            %TOTAL Calculate total mass balance
            % spec = concentration matrix (time x sections)
            % operators = struct with coagulation and linear operators
            % Returns: total_gains, total_losses structs

            [n_times, n_sections] = size(spec);

            % Total sinking losses
            sinking = diag(operators.sink_loss);
            sinking = sinking';
            sinking = sinking(ones(n_times, 1), :);
            sink_losses = sinking .* spec;
            total_sink_losses = sum(sink_losses, 2);

            % Net growth
            net_growth = zeros(n_times, 1);
            for i_time = 1:n_times
                v = spec(i_time, :)';
                g1 = operators.growth * v;
                net_growth(i_time) = sum(g1);
            end

            % Coagulation losses from system (largest section)
            coag_losses = zeros(n_times, 1);
            for i_time = 1:n_times
                v = spec(i_time, :);
                v_r = v';

                % Loss from largest section
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

            % Assemble results
            total_gains.growth = net_growth;
            total_losses.sett = total_sink_losses;
            total_losses.coag = coag_losses;
        end

        function [gains, losses] = sectionalWithRates(spec, rhs)
            %SECTIONALWITHRATES Calculate sectional mass balance using RHS
            % Alternative method using CoagulationRHS object
            % spec = concentration matrix (time x sections)
            % rhs = CoagulationRHS object
            % Returns: gains, losses structs

            [n_times, n_sections] = size(spec);

            % Initialize arrays
            coag_gains = zeros(n_times, n_sections);
            coag_losses = zeros(n_times, n_sections);

            % Calculate using rate terms from RHS
            for i_time = 1:n_times
                [term1, term2, term3, term4, term5] = rhs.rateTerms(spec(i_time, :)');

                coag_gains(i_time, :) = term1 + term2;
                coag_losses(i_time, :) = -term4;  % Disaggregation loss
            end

            % For now, return basic structure
            % TODO: Implement full sectional analysis with RHS
            gains.coag = coag_gains;
            gains.growth = zeros(n_times, n_sections);  % Placeholder

            losses.coag = coag_losses;
            losses.growth = zeros(n_times, n_sections);  % Placeholder
            losses.settl = zeros(n_times, n_sections);   % Placeholder
        end

        function displayBalanceSummary(gains, losses, time_vector)
            %DISPLAYBALANCESUMMARY Display summary of mass balance
            % gains, losses = structs from sectional or total analysis
            % time_vector = time points

            fprintf('Mass Balance Summary:\n');
            fprintf('  Time points: %d\n', length(time_vector));
            fprintf('  Time range: [%.2f, %.2f]\n', time_vector(1), time_vector(end));

            % Display gains
            if isfield(gains, 'coag')
                fprintf('  Coagulation gains range: [%.2e, %.2e]\n', ...
                    min(gains.coag(:)), max(gains.coag(:)));
            end
            if isfield(gains, 'growth')
                fprintf('  Growth gains range: [%.2e, %.2e]\n', ...
                    min(gains.growth(:)), max(gains.growth(:)));
            end

            % Display losses
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
