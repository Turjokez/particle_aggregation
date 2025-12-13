function simulation_runner_3
% SIMULATION_RUNNER_3
%
% Adrian mass-balance tests for the coagulation model.
%
% TEST 1: Constant forcing, PP ON
%   - Keep epsilon, mu, and PP constant.
%   - Run to steady state.
%   - Compute:
%       PP(t) = mu * sum_z C1(z,t) * V1 * dz
%     where C1 is the smallest size class.
%   - At steady state, check whether:
%       <PP>_SS ≈ <F_bottom>_SS
%     for three values of mu: 0.1, 0.5, 1.0 d^-1.
%
% TEST 2: PP OFF, no other biology
%   - Switch off PP.
%   - Keep turbulence and everything else constant.
%   - Initial condition: particles only in the top layer.
%   - Run 30 days.
%   - Check mass balance:
%       M0 ≈ M_final + ∫ F_bottom dt
%     where M is the total biovolume in 0–200 m.
%
% Units:
%   - Concentrations in result.concentrations are in number / volume.
%   - We multiply by particle volume V_m3 to get biovolume, and by dz
%     to integrate over depth, so all diagnostics are consistent.

clc;
fprintf('=============================================\n');
fprintf('   COAGULATION MODEL – MASS BALANCE TESTS\n');
fprintf('=============================================\n\n');

%% COMMON SETTINGS
z_max = 200;
dz    = 10;
eps0  = 1e-7;  % constant turbulence (W kg^-1)

outDir = 'fig_massTests';
if ~exist(outDir,'dir'); mkdir(outDir); end

%% TEST 1: Constant forcing with PP ON (steady state check)
fprintf('TEST 1: Constant forcing with PP ON (steady state check)\n');

mu_vals = [0.10, 0.50, 1.00];  % d^-1

for ii = 1:numel(mu_vals)
    mu_val = mu_vals(ii);
    fprintf('  -- Running mu = %.2f d^-1\n', mu_val);

    T_final = 60;                % days
    dt      = 0.1;
    t_run   = (0:dt:T_final).';

    % ---------- CONFIG ----------
    cfg = SimulationConfig();
    cfg.use_column = true;
    cfg.z_max      = z_max;
    cfg.dz         = dz;

    % Turn off extra biology so we are mainly looking at
    % NPP + attenuation + coagulation + sinking.
    cfg.growth      = 0.0;
    cfg.gro_sec     = 0;
    cfg.loss_rate   = 0.0;
    cfg.disagg_use_nonlinear = false;

    % Attenuation: depth-independent mu
    cfg.attenuation_rate         = mu_val;   % mu (d^-1)
    cfg.attenuation_depth_factor = 0.0;

    % NPP: keep as whatever your model uses; we only diagnose PP
    % using Adrian's definition (mu * sum C1).
    cfg.use_NPP         = true;
    cfg.NPP_rate        = mu_val;          % simple mapping
    cfg.NPP_t_step      = Inf;             % always on
    cfg.NPP_rate_after  = cfg.NPP_rate;

    % Time-stepping options inside config
    cfg.t_init  = 0;
    cfg.t_final = T_final;
    cfg.delta_t = dt;

    % Constant epsilon time series
    eps_run = eps0 * ones(size(t_run));

    % ---------- RUN SIMULATION ----------
    sim = CoagulationSimulation(cfg);
    sim.setEpsilonTimeSeries(t_run, eps_run);
    sim.run('tspan', t_run);   % let the class choose its default v0

    result   = sim.result;
    grid_obj = sim.grid;

    % ---------- DIAGNOSTICS ----------
    [F_bottom, PP_in, M_total] = ...
        diagnostics_test1(result, grid_obj, cfg);

    % steady-state averages (last 10% of time)
    Nt     = numel(result.time);
    idx_ss = round(0.9*Nt):Nt;

    F_ss  = mean(F_bottom(idx_ss));
    PP_ss = mean(PP_in(idx_ss));

    ratio = F_ss / PP_ss;

    fprintf('Simulation Complete.\n');
    fprintf('     <F_bottom>_SS = %.4e (mass / day)\n', F_ss);
    fprintf('     <PP input>_SS = %.4e (mass / day)\n', PP_ss);
    fprintf('     Ratio F_out / PP = %.3f\n\n', ratio);

    % ---------- PLOT ----------
    fig = figure('Color','w','Position',[80 80 1200 700]);
    t = result.time(:);

    % 1) Flux out
    subplot(3,1,1);
    plot(t, F_bottom, 'b-', 'LineWidth', 1.5);
    ylabel('Flux out');
    title(sprintf('Test 1,  \\mu = %.2f  d^{-1}', mu_val));
    grid on;
    xlim([t(1) t(end)]);

    % 2) PP input
    subplot(3,1,2);
    plot(t, PP_in, 'r-', 'LineWidth', 1.5);
    ylabel('PP input');
    title('Total PP input (biovolume / area / day)');
    grid on;
    xlim([t(1) t(end)]);

    % 3) Mass in column
    subplot(3,1,3);
    plot(t, M_total, 'k-', 'LineWidth', 1.5);
    ylabel('Mass in column');
    xlabel('Time (days)');
    title('Total mass (biovolume) in 0–200 m');
    grid on;
    xlim([t(1) t(end)]);

    fname_base = sprintf('Test1_mu_%0.2f', mu_val);
    saveas(fig, fullfile(outDir, [fname_base '.png']));
    savefig(fig, fullfile(outDir, [fname_base '.fig']));
end

fprintf('TEST 1 complete.\n\n');

%% TEST 2: PP OFF, NO BIOLOGY (pure coag + sinking)
fprintf('TEST 2: PP OFF, no biology (mass-conservation test)\n');

T_final = 30;
dt      = 0.1;
t_run   = (0:dt:T_final).';

cfg2 = SimulationConfig();
cfg2.use_column = true;
cfg2.z_max      = z_max;
cfg2.dz         = dz;

% Turn OFF all explicit sources/sinks other than coag + sinking
cfg2.use_NPP         = false;
cfg2.NPP_rate        = 0;
cfg2.NPP_t_step      = 0;
cfg2.NPP_rate_after  = 0;

cfg2.attenuation_rate         = 0.0;
cfg2.attenuation_depth_factor = 0.0;

cfg2.growth      = 0.0;
cfg2.gro_sec     = 0;
cfg2.loss_rate   = 0.0;

cfg2.disagg_use_nonlinear = false;

% time range in config (for internal defaults)
cfg2.t_init  = 0;
cfg2.t_final = T_final;
cfg2.delta_t = dt;

% constant epsilon
eps_run2 = eps0 * ones(size(t_run));

sim2 = CoagulationSimulation(cfg2);
sim2.setEpsilonTimeSeries(t_run, eps_run2);

% --- Initial condition: particles only in top layer, smallest bin ---
Nz = floor(cfg2.z_max / cfg2.dz);
Ns = cfg2.n_sections;

Y0 = zeros(Nz*Ns, 1);
C0 = 1e-6;            % base concentration in smallest bin
Y0(1) = C0;           % top layer, smallest size class

% Run with explicit v0
sim2.run('tspan', t_run, 'v0', Y0);

result2   = sim2.result;
grid_obj2 = sim2.grid;

[F_bottom_2, M_total_2, M0, Mexport_int, mass_balance] = ...
    diagnostics_test2(result2, grid_obj2, cfg2, t_run);

fprintf('Simulation Complete.\n');
fprintf('  Initial mass in column      = %.4e\n', M0);
fprintf('  Final mass in column        = %.4e\n', M_total_2(end));
fprintf('  Integrated exported mass    = %.4e\n', Mexport_int);
fprintf('  Final + Exported - Initial  = %.4e (should be ~ 0)\n\n', ...
        mass_balance);

fig2 = figure('Color','w','Position',[100 100 1000 600]);
t2 = result2.time(:);

subplot(2,1,1);
plot(t2, M_total_2, 'k-','LineWidth',1.5);
ylabel('Mass in column');
title('Test 2: Mass vs Time');
grid on;
xlim([t2(1) t2(end)]);

subplot(2,1,2);
plot(t2, F_bottom_2, 'b-','LineWidth',1.5);
ylabel('Flux out bottom');
xlabel('Time (days)');
title('Test 2: Flux out bottom vs Time');
grid on;
xlim([t2(1) t2(end)]);

fname2 = 'Test2_mass_flux';
saveas(fig2, fullfile(outDir, [fname2 '.png']));
savefig(fig2, fullfile(outDir, [fname2 '.fig']));

fprintf('TEST 2 complete.\nDone.\n');

end

%% ---- Helper for TEST 1 (Adrian PP definition) ----
function [F_bottom, PP_input, M_total] = diagnostics_test1(result, grid_obj, cfg)

t  = result.time(:);
Y  = result.concentrations;    % [Nt x (Nz*Ns)]
Nt = numel(t);
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;

% particle geometry
r_cm    = grid_obj.getFractalRadii();
r_v     = grid_obj.getConservedRadii();
V_m3    = (4/3)*pi*r_v.^3 * 1e-6;   % cm^3 -> m^3
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_obj.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;  % cm/s -> m/day

dz = cfg.dz;

F_bottom = zeros(Nt,1);
PP_input = zeros(Nt,1);
M_total  = zeros(Nt,1);

mu = cfg.attenuation_rate;  % μ in d^-1

for it = 1:Nt
    Yt = Y(it,:);                     % 1 x (Nz*Ns)
    Yz = reshape(Yt, [Ns, Nz]).';     % [Nz x Ns]

    % ---- Flux at bottom (200 m) ----
    Yb = Yz(Nz,:);  % bottom layer concentrations
    % biovolume flux out = sum( C * V * w )
    F_bottom(it) = sum( Yb .* V_m3.' .* W_m_d.' );   % mass/area/day

    % ---- PP input (Adrian: μ * sum_z C1(z,t)) ----
    C1 = Yz(:,1);   % smallest size class across depth
    PP_input(it) = mu * sum( C1 .* V_m3(1) * dz );

    % ---- Total mass in column ----
    % biovolume per area: sum_z sum_s C(z,s) * V_s * dz
    M_total(it) = sum( Yz .* V_m3.', 'all' ) * dz;
end

end

%% ---- Helper for TEST 2 ----
function [F_bottom, M_total, M0, Mexport_int, mass_balance] = ...
          diagnostics_test2(result, grid_obj, cfg, t_run)

t  = result.time(:);
Y  = result.concentrations;
Nt = numel(t);
Nz = floor(cfg.z_max / cfg.dz);
Ns = cfg.n_sections;

r_cm    = grid_obj.getFractalRadii();
r_v     = grid_obj.getConservedRadii();
V_m3    = (4/3)*pi*r_v.^3 * 1e-6;
ws_cm_s = SettlingVelocityService.velocity(r_cm, r_v, grid_obj.setcon);
W_m_d   = (ws_cm_s / 100) * 86400;

dz = cfg.dz;

F_bottom = zeros(Nt,1);
M_total  = zeros(Nt,1);

for it = 1:Nt
    Yt = Y(it,:);
    Yz = reshape(Yt, [Ns, Nz]).';

    Yb = Yz(Nz,:);
    F_bottom(it) = sum( Yb .* V_m3.' .* W_m_d.' );

    M_total(it) = sum( Yz .* V_m3.', 'all' ) * dz;
end

% initial mass from first time step
M0 = M_total(1);

% integrate exported mass over time (simple trapezoid)
Mexport_int = trapz(t_run, F_bottom);   % same units as M_total

mass_balance = M_total(end) + Mexport_int - M0;

end