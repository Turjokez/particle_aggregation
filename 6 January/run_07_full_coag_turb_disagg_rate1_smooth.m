% run_07_full_coag_turb_disagg_rate1_smooth.m
clear; close all;

cfg = SimulationConfig();
cfg.use_column  = true;
cfg.z_max       = 65;
cfg.dz          = 5;
cfg.n_sections  = 24;

cfg.enable_coag = true;

if isprop(cfg,'enable_disagg');     cfg.enable_disagg = true; end
if isprop(cfg,'disagg_apply_in');   cfg.disagg_apply_in = 'rhs'; end

% disagg rate (start small)
if isprop(cfg,'disagg_rate')
    cfg.disagg_rate = 1.0;   % [d^-1]
end

cfg.growth_mode = 'shift';

% turbulence pulse (SMOOTH)
eps0 = 1e-8;
eps1 = 1e-5;
t1   = 18;
t2   = 20;
tau  = 0.05; % days smoothing

if isprop(cfg,'eps_fun')
    cfg.eps_fun = @(t,z) (eps0 + (eps1-eps0) .* ...
        (0.5*(tanh((t - t1)/tau) - tanh((t - t2)/tau))));
end

if isprop(cfg,'eps_ref')
    cfg.eps_ref = eps0;
end

S = CoagulationSimulation(cfg);
out = S.run();

% STOP if solver failed early
if out.time(end) < cfg.t_final - 1e-6
    error('Solver did not reach t_final. No diagnostics saved.');
end

DiagnosticsSuite.runPart1(S, 'run_07_rate1_smooth');
DiagnosticsSuite.runPart2(S, 'run_07_rate1_smooth');

disp('RUN 07 DONE');