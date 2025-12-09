%% run_with_obs_epsilon.m
% Driver: load Adrian's ε(t), run the model, make all figures.

clear; close all; clc; clear classes; rehash;  % refresh class defs

% --- 0) (optional) add your repo to the path
% addpath('/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code');

% --- 1) Pick Adrian's file (or hard-code the path)
matfile = '/Users/turjo/Library/CloudStorage/OneDrive-UniversityofGeorgia/particle_aggregation/test_code/epsilon_daily.mat';
% [fn,fp] = uigetfile('*.mat','Select epsilon_daily.mat', fileparts(matfile));
% if isequal(fn,0), error('No file selected.'); end
% matfile = fullfile(fp,fn);

% --- 2) Load ML-mean ε(t) (linear units, daily, light smoothing)
[t_days, eps_series] = load_epsilon_series(matfile, 0, ...
    'time_field','mtime', ...
    'eps_field','eps', ...
    'depth_field','z', ...
    'mld_field','mld', ...
    'unit','linear', ...          % IMPORTANT: your file is already linear
    'agg','mld', 'topN',3, ...
    'resample_daily',true, 'smooth_win',3);

fprintf('ε range: %.2e to %.2e W kg^-1\n', min(eps_series), max(eps_series));

% --- 3) Configure simulation (keep all your physics the same)
cfg = SimulationConfig( ...
    'epsilon_profile','observed', ...   % << use the observed series
    'epsilon_time',   t_days, ...
    'epsilon_series', eps_series, ...
    'epsilon_ref',    1e-6, ...
    'disagg_use_nonlinear', true);      % same as your working runs

% --- 4) Run and plot (OutputGenerator now has guardrails)
sim = CoagulationSimulation(cfg);
result = sim.run();
sim.generateOutputs(true);              % Figures 1–4