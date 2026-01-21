% ==============================================================
% FILE: run_06_coag_only_0D.m
% True coag-only conservation test (no solver hacks)
% ==============================================================
clear; clc;

cfg = SimulationConfig();

% --------------------------------------------------------------
% 0-D
% --------------------------------------------------------------
cfg.use_column = false;

% --------------------------------------------------------------
% Physics switches
% --------------------------------------------------------------
cfg.enable_coag    = true;
cfg.enable_disagg  = false;
cfg.enable_pp      = false;

% No growth / no linear / no sinking
cfg.growth_mode    = 'shift';
cfg.growth         = 0;
cfg.gro_sec        = 0;

cfg.enable_linear  = false;
cfg.enable_sinking = false;

% --------------------------------------------------------------
% State convention
% --------------------------------------------------------------
cfg.state_is_biovolume      = false;   % NUMBER state
cfg.enforce_coag_bv_closure = false;   % IMPORTANT: OFF for operator test

% --------------------------------------------------------------
% Beta orientation fix (your NEW switch)
% --------------------------------------------------------------
cfg.beta_row_fix           = true;     % must be ON to activate fix
cfg.beta_transpose_b1b3    = true;     % transpose b1 and b3 only

% --------------------------------------------------------------
% IMPORTANT: turn off safety hacks for a pure operator test
% --------------------------------------------------------------
cfg.use_nonnegative = false;   % do not let solver enforce positivity
cfg.clip_negative   = false;   % do not clip inside RHS

% --------------------------------------------------------------
% Run
% --------------------------------------------------------------
sim = CoagulationSimulation(cfg);
out = sim.run();

t  = out.time(:);
Y  = out.concentrations;     % Nt x Ns (NUMBER)
av = sim.grid.av_vol(:);     % Ns x 1

% --------------------------------------------------------------
% Two "mass" checks (keep BOTH)
%   1) Physical biovolume inventory: sum(n * av)
%   2) Operator-null weight found by brute force: sum(n * av.^2)
% --------------------------------------------------------------
w_phys = av;
w_op   = av.^2;

M_phys = Y * w_phys;
M_op   = Y * w_op;

fprintf('\n=== COAG-ONLY 0-D CONSERVATION TEST ===\n');

fprintf('[PHYS] M(t0)=%.6e  M(tend)=%.6e  rel=%.3e   (biovolume-weight)\n', ...
    M_phys(1), M_phys(end), (M_phys(end)-M_phys(1))/max(M_phys(1),realmin));

fprintf('[OP  ] M(t0)=%.6e  M(tend)=%.6e  rel=%.3e   (av_vol.^2 weight)\n', ...
    M_op(1), M_op(end), (M_op(end)-M_op(1))/max(M_op(1),realmin));

% --------------------------------------------------------------
% Instantaneous leak check
% --------------------------------------------------------------
it = max(2, round(numel(t)*0.6));
v  = Y(it,:).';
dv = sim.rhs.rhs(t(it), v);

fprintf('\nCHECK at t=%.3f d\n', t(it));
fprintf('||dv||_inf               = %.6e\n', norm(dv,inf));
fprintf('[PHYS] sum(dv .* av_vol) = %.6e  (should be ~0 if biovolume conserved)\n', sum(dv .* w_phys));
fprintf('[OP  ] sum(dv .* av^2)   = %.6e  (should be ~0 if av^2 is conserved)\n', sum(dv .* w_op));
fprintf('======================================\n');