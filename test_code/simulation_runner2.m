function simulation_runner2
% SIMULATION_RUNNER2
%
% Adrian mass-balance tests (no EOFs):
%   Test 1: Constant forcing, with primary production (PP) on.
%           Run to steady state and compare:
%               <flux out of bottom>  ~  <PP input>
%
%   Test 2: PP off, start with particles only in the surface layer.
%           Run ~30 days and check:
%               initial mass  ~  exported mass + final mass
%
% This script does NOT change any of your class files.

clc;
fprintf('=============================================\n');
fprintf('   COAGULATION MODEL – MASS BALANCE TESTS\n');
fprintf('=============================================\n\n');

%% ----- USER SETTINGS (same for both tests) -------------------------
Nz          = 20;       % 200 m / 10 m
zMax        = 200;
dz          = 10;
eps0        = 1e-7;     % constant epsilon (W kg^-1)
mu_list     = [0.1 0.5 1.0];   % d^-1, used as PP rate in Test 1
T_steady    = 60;       % days for Test 1
T_decay     = 30;       % days for Test 2
dt          = 0.1;      % days

%% ======================= TEST 1 ====================================
fprintf('TEST 1: Constant forcing with PP ON (steady state check)\n');
fprintf('        eps = %.1e W kg^-1, run %.1f days\n\n', eps0, T_steady);

for imu = 1:numel(mu_list)
    mu = mu_list(imu);

    fprintf('  -- Running mu = %.2f d^-1\n', mu);

    % --- build config ---
    cfg = SimulationConfig();
    cfg.use_column = true;
    cfg.z_max      = zMax;
    cfg.dz         = dz;

    % no attenuation here – just PP and sinking
    cfg.attenuation_rate         = 0.0;
    cfg.attenuation_depth_factor = 0.0;

    % use NPP as simple PP in the first size class
    cfg.use_NPP        = true;
    cfg.NPP_rate       = mu;           % acts like mu * C1
    cfg.NPP_profile    = 'constant';
    cfg.NPP_section    = 1;
    cfg.NPP_t_step     = 0;            % no change in time
    cfg.NPP_rate_after = mu;

    % time axis and constant epsilon
    t_run   = (0:dt:T_steady).';
    eps_run = eps0 * ones(size(t_run));

    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    sim.run('tspan', t_run);

    % --- diagnostics: flux and PP input ---
    [F_bottom, mass_col, NPP_in] = ...
        diagnostics_test1(sim.result, sim.grid, cfg, mu);

    % pick "steady state" as the last third of the run
    Nt     = numel(t_run);
    i_ss   = round(2*Nt/3) : Nt;
    F_ss   = mean(F_bottom(i_ss));
    PP_ss  = mean(NPP_in(i_ss));

    fprintf('     <F_bottom>_SS = %.4e (mass / day)\n', F_ss);
    fprintf('     <PP input>_SS = %.4e (mass / day)\n', PP_ss);
    fprintf('     Ratio F_out / PP = %.3f\n\n', F_ss / PP_ss);

    % small sanity plot (optional)
    figure('Name',sprintf('Test1_mu%.2f',mu),'Color','w');
    subplot(3,1,1);
    plot(t_run, F_bottom, 'b-','LineWidth',1.5);
    ylabel('Flux out');
    title(sprintf('Test 1, \\mu = %.2f d^{-1}',mu));
    grid on;

    subplot(3,1,2);
    plot(t_run, NPP_in, 'r-','LineWidth',1.5);
    ylabel('PP input');
    grid on;

    subplot(3,1,3);
    plot(t_run, mass_col, 'k-','LineWidth',1.5);
    ylabel('Mass in column');
    xlabel('Time (days)');
    grid on;
end

fprintf('TEST 1 complete.\n\n');

%% ======================= TEST 2 ====================================
fprintf('TEST 2: PP OFF, initial particles only in top layer\n');
fprintf('        eps = %.1e W kg^-1, run %.1f days\n\n', eps0, T_decay);

% --- build config ---
cfg2 = SimulationConfig();
cfg2.use_column = true;
cfg2.z_max      = zMax;
cfg2.dz         = dz;

cfg2.use_NPP        = false;   % PP off
cfg2.attenuation_rate         = 0.0;
cfg2.attenuation_depth_factor = 0.0;

% time and epsilon
t_run2   = (0:dt:T_decay).';
eps_run2 = eps0 * ones(size(t_run2));

% initial condition: particles only in top layer
grid2    = cfg2.derive();
v0_slab  = InitialSpectrumBuilder.initialSpectrum(cfg2, grid2);
Ns       = cfg2.n_sections;
Nz       = floor(cfg2.z_max / cfg2.dz);

v0_grid       = zeros(Nz, Ns);
v0_grid(1,:)  = v0_slab;   % only top layer
v0_vec        = v0_grid.'; % [Ns x Nz]
v0_vec        = v0_vec(:); % [Nz*Ns x 1]

sim2 = CoagulationSimulation(cfg2);
sim2.setEpsilonTimeSeries(t_run2, eps_run2);
sim2.run('tspan', t_run2, 'v0', v0_vec);

[F_bottom2, mass_col2] = diagnostics_test2(sim2.result, sim2.grid, cfg2);

% mass balance numbers
initial_mass  = mass_col2(1);
final_mass    = mass_col2(end);
exported_mass = trapz(t_run2, F_bottom2);   % integral of flux over time

fprintf('  Initial mass in column      = %.4e\n', initial_mass);
fprintf('  Final mass in column        = %.4e\n', final_mass);
fprintf('  Integrated exported mass    = %.4e\n', exported_mass);
fprintf('  Final + Exported - Initial  = %.4e (should be ~ 0)\n\n', ...
        final_mass + exported_mass - initial_mass);

% simple plot
figure('Name','Test2_decay','Color','w');
subplot(2,1,1);
plot(t_run2, mass_col2, 'k-','LineWidth',1.5);
ylabel('Mass in column');
title('Test 2: Mass vs Time');
grid on;

subplot(2,1,2);
plot(t_run2, F_bottom2, 'b-','LineWidth',1.5);
ylabel('Flux out bottom');
xlabel('Time (days)');
grid on;

fprintf('TEST 2 complete.\n');
fprintf('Done.\n');

end

%% ================== DIAGNOSTIC HELPERS ===========================

function [F_bottom, mass_col, NPP_in] = diagnostics_test1(result, grid, cfg, mu)
% Bottom flux, column mass, and simple PP input estimate for Test 1.

t  = result.time(:); %#ok<NASGU> % not used directly here, but kept for clarity
Y  = result.concentrations;
Nt = size(Y,1);
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;

% radii and velocities
r_cm    = grid.getFractalRadii();
r_v     = grid.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;              % cm/s -> m/day
V_m3    = (4/3) * pi * r_v.^3 * 1e-6;           % cm^3 -> m^3

% bottom flux
export_depth_idx = Nz;
Q_bottom = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);
F_bottom = sum(Q_bottom .* V_m3.' .* W_m_d.', 2);

% column mass at each time
mass_col = zeros(Nt,1);
for k = 1:Nz
    idx = (k-1)*Ns + (1:Ns);
    Q_layer = Y(:,idx);
    mass_col = mass_col + sum(Q_layer .* V_m3.', 2);
end

% PP input estimate: mu * sum(C1 over depth)
idx_first = 1:Ns:(Nz*Ns);  % first size class in each layer
C1_tot    = sum(Y(:, idx_first), 2);
NPP_in    = mu * C1_tot;   % proportional to total C1

end


function [F_bottom, mass_col] = diagnostics_test2(result, grid, cfg)
% Bottom flux and column mass for Test 2 (no PP).

Y  = result.concentrations;
Nt = size(Y,1);
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;

% radii and velocities
r_cm    = grid.getFractalRadii();
r_v     = grid.getConservedRadii();
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;              % cm/s -> m/day
V_m3    = (4/3) * pi * r_v.^3 * 1e-6;           % cm^3 -> m^3

% bottom flux
export_depth_idx = Nz;
Q_bottom = Y(:, (export_depth_idx-1)*Ns + 1 : export_depth_idx*Ns);
F_bottom = sum(Q_bottom .* V_m3.' .* W_m_d.', 2);

% total mass in column
mass_col = zeros(Nt,1);
for k = 1:Nz
    idx = (k-1)*Ns + (1:Ns);
    Q_layer = Y(:,idx);
    mass_col = mass_col + sum(Q_layer .* V_m3.', 2);
end

end