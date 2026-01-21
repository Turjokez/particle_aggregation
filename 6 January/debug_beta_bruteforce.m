% ==============================================================
% FILE: debug_beta_bruteforce.m
% Purpose: brute-force beta transpose flags at t0 and report best
% Requires: you already ran sim = CoagulationSimulation(cfg); out=sim.run();
% ==============================================================

% --- sanity: beta flag really on ---
fprintf('cfg.beta_row_fix (cfg)        = %d\n', cfg.beta_row_fix);
fprintf('cfg.beta_row_fix (rhs.config) = %d\n', sim.rhs.config.beta_row_fix);

Ns = cfg.n_sections;

Y  = out.concentrations;       % Nt x Ns
n0 = Y(1,:).';                 % NUMBER state at t0 (Ns x 1)
nr = n0.';                     % row (1 x Ns)
nshift = [0, nr(1:Ns-1)];      % shifted row
av = sim.grid.av_vol(:);       % Ns x 1

B = sim.rhs.betas;

names = {'b1','b2','b3','b4','b5'};
Ms = {B.b1, B.b2, B.b3, B.b4, B.b5};

best_abs = Inf;
best_flags = [];
best_leak = NaN;

fprintf('\n=== BRUTE FORCE BETA ORIENTATION CHECK (t0) ===\n');

for mask = 0:(2^5-1)
    tr = bitget(mask,1:5); % [b1 b2 b3 b4 b5] 1=transpose
    leak = test_one(tr, nr, nshift, av, Ms);
    absleak = abs(leak);

    if absleak < best_abs
        best_abs   = absleak;
        best_flags = tr;
        best_leak  = leak;
    end
end

fprintf('BEST abs(leak) = %.6e  leak=%.6e\n', best_abs, best_leak);
fprintf('Transpose flags [b1 b2 b3 b4 b5] = [%d %d %d %d %d]\n', best_flags);

fprintf('Meaning:\n');
for i=1:5
    if best_flags(i)==1
        fprintf('  %s: use nr*(%s'')\n', names{i}, names{i});
    else
        fprintf('  %s: use nr*%s\n', names{i}, names{i});
    end
end
fprintf('==============================================\n');

% ==============================================================
% Local function MUST be at end of file (MATLAB rule)
% ==============================================================
function leak = test_one(tr, nr, nshift, av, Ms)
    % tr(i)=1 means transpose Mi
    out = cell(1,5);
    for ii=1:5
        Mi = Ms{ii};
        if tr(ii)==1
            out{ii} = nr*(Mi.');
        else
            out{ii} = nr*Mi;
        end
    end

    % same algebra as RHS
    tg2 = nr .* out{2};
    tl3 = nr .* out{3};
    tl4 = nr .* out{4};
    tl5 = nr .* out{5};

    term1 = tg2 - tl3 - tl4 - tl5;
    term2 = out{1} .* nshift;

    dn = (term1 + term2).';  % Ns x 1
    leak = sum(dn .* av);
end