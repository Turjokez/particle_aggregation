% RUN_PART1_DIAGNOSTICS

cfg = SimulationConfig();
cfg.use_column = true;
cfg.z_max = 65;
cfg.dz = 5;

cfg.n_sections = 24;   % so we can have D >= 2000 um

% ==========================================================
% NEW-2026-01-07: one-switch mode control
%   debug_linear_only = true  -> sinking + linear only (no coag, no disagg)
%   debug_linear_only = false -> full physics (coag ON; disagg optional)
% ==========================================================
debug_linear_only = false;   % <<< change this only

% ----------------------------------------------------------
% Coagulation ON/OFF (safe even if property doesn't exist)
% ----------------------------------------------------------
if isprop(cfg,'enable_coag')
    if debug_linear_only
        cfg.enable_coag = false;
    else
        cfg.enable_coag = true;
    end
end

% ----------------------------------------------------------
% Disaggregation / turbulence forcing (optional)
%   Keep OFF for now unless you want turbulence experiments.
% ----------------------------------------------------------
if isprop(cfg,'enable_disagg')
    if debug_linear_only
        cfg.enable_disagg = false;
    else
        cfg.enable_disagg = false;   % set true if you want breakup physics
    end
end

% If you DO enable disagg, you need epsilon provided:
% Option 1: constant epsilon
% if isprop(cfg,'epsilon_const')
%     cfg.epsilon_const = 1e-7;   % example
% end
%
% Option 2: time series epsilon (surface or Nt x Nz)
% cfg.epsilon_time   = [...]; 
% cfg.epsilon_series = [...];

% ----------------------------------------------------------
% Run
% ----------------------------------------------------------
S = CoagulationSimulation(cfg);
out = S.run();

DiagnosticsSuite.runPart1(S);  % saves to ./diagnostics_part1/