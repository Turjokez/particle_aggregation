% quick_run_epsilon.m
clear; clc;

addpath(genpath('/path/to/your/project/root'));  % <-- update

% --- load observed ε (same helper you used) ---
matfile = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code/epsilon_daily.mat';
[t_days, eps_series] = load_epsilon_series(matfile,0, ...
    'time_field','mtime','eps_field','eps','depth_field','z','mld_field','mld', ...
    'unit','linear','agg','mld','topN',3,'resample_daily',true,'smooth_win',3);

% IMPORTANT: floor ε to avoid zeros (prevents α → ∞ when we enable scaling)
eps_floor   = 1e-8;                     % W kg^-1
eps_series  = max(eps_series, eps_floor);

% --- config: observed forcing, but NO α(ε) scaling yet ---
cfg = SimulationConfig( ...
    'epsilon_profile','observed', ...
    'epsilon_time',t_days, ...
    'epsilon_series',eps_series, ...
    'epsilon_ref',1e-6, ...
    'p_alpha',0.0, ...                 % turn OFF stickiness scaling for this check
    'disagg_use_nonlinear', true, ...
    't_final',30, 'delta_t',1);

sim = CoagulationSimulation(cfg);
res = sim.run();            % integrates
sim.generateOutputs(true);  % plots + console summary

% quick check on ε actually used
fprintf('ε used: min=%.2e, median=%.2e, max=%.2e\n', ...
    min(res.epsilon_used), median(res.epsilon_used), max(res.epsilon_used));