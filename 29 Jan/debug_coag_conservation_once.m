% debug_coag_conservation_once.m
clear; clc;

cfg = SimulationConfig();

% Force 0-D
if isprop(cfg,'use_column'), cfg.use_column = false; end

% Turn ON only coag
if isprop(cfg,'enable_coag'),    cfg.enable_coag = true;  end
if isprop(cfg,'enable_disagg'),  cfg.enable_disagg = false; end
if isprop(cfg,'enable_sinking'), cfg.enable_sinking = false; end
if isprop(cfg,'enable_pp'),      cfg.enable_pp = false; end
if isprop(cfg,'growth'),         cfg.growth = 0; end

% Make it super fast (just to initialize internal objects)
if isprop(cfg,'t_max'),   cfg.t_max = 0.01; end
if isprop(cfg,'dt_out'),  cfg.dt_out = 0.01; end

sim = CoagulationSimulation(cfg);

% ---- KEY: run once so rhs gets built in your code version ----
out = sim.run(); %#ok<NASGU>

% Now rhs should exist
rhs = sim.rhs;

fprintf('\nDEBUG: class(sim.rhs) = %s\n', class(rhs));
if isempty(rhs)
    error('sim.rhs is still empty after sim.run(). Your version stores RHS differently.');
end

% Try to get section grid
if ismethod(rhs,'getSectionGrid')
    gridS = rhs.getSectionGrid();
else
    error('rhs does not have getSectionGrid(). We need to locate where the section grid lives in your version.');
end

Ns = cfg.n_sections;

% Positive random state
v = abs(rand(Ns,1)) + 1e-12;

% Bin volumes weight
Vbin = Disaggregation.getBinVolumes_cm3(gridS);

% Coag tendency
T = rhs.decomposeTerms(0.0, v);
dv_coag = T.dv_coag(:);

leak_sum = sum(dv_coag);
leak_V   = sum(dv_coag .* Vbin);

scale = max(1e-30, norm(dv_coag,1));

fprintf('\n--- COAG CONSERVATION CHECK (single state) ---\n');
fprintf('sum(dv_coag)                 = %.3e\n', leak_sum);
fprintf('sum(dv_coag .* Vbin)         = %.3e\n', leak_V);
fprintf('relative leak (Vbin)         = %.3e\n', abs(leak_V)/scale);
fprintf('--------------------------------------------\n');