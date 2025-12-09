clear; clc;
cfg = SimulationConfig( ...
   'epsilon_profile','constant', ...   % no time variation
   'epsilon',1e-7, ...
   'disagg_use_nonlinear',false, ...   % external breakup OFF
   'p_alpha',0, ...
   'use_NPP',false ...
   );
% If your old baseline did not include growth, uncomment the next line:
% cfg.growth = 0;

sim = CoagulationSimulation(cfg);
res = sim.run();
sim.generateOutputs(true);