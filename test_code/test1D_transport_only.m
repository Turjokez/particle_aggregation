function test1D_transport_only()
% TEST1D_TRANSPORT_ONLY
% 1-D column: sinking transport + coagulation ONLY.
% Checks column mass budget:  Mf + E - M0 ≈ 0

clc;
fprintf('=============================================\n');
fprintf('   TEST 1-D: TRANSPORT + COAG (NO BIOLOGY)\n');
fprintf('=============================================\n');

%% --- CONFIG ----------------------------------------------------------
cfg = SimulationConfig();
cfg.use_column   = true;
cfg.n_sections   = 35;
cfg.z_max        = 200;          % [m]
cfg.dz           = 10;           % [m]

cfg.t_init       = 0.0;
cfg.t_final      = 30.0;         % [days]
cfg.delta_t      = 0.1;

% physics (mild coag)
cfg.epsilon_profile = 'constant';
cfg.epsilon         = 1e-9;
cfg.alpha_base      = 0.01;

% *** TURN OFF ALL BIOLOGY / LOSSES ***
cfg.use_NPP                  = false;
cfg.growth                   = 0.0;
cfg.attenuation_rate         = 0.0;
cfg.attenuation_depth_factor = 0.0;
cfg.loss_rate                = 0.0;

% no nonlinear breakup
cfg.disagg_use_nonlinear     = false;

%% --- GRID & INITIAL COLUMN ------------------------------------------
dg   = cfg.derive();
Ns   = cfg.n_sections;
Nz   = floor(cfg.z_max / cfg.dz);

v_vol = dg.av_vol(:);            % [cm^3]
dz_cm = cfg.dz * 100;            % [cm]

% slab spectrum (used only for top layer)
N0_slab = InitialSpectrumBuilder.initialSpectrum(cfg, dg);

N0_grid      = zeros(Nz, Ns);
N0_grid(1,:) = N0_slab;          % only top layer filled

Y0 = N0_grid.';                  % [Ns x Nz]
Y0 = Y0(:);                      % [Ns*Nz x 1]

% initial column mass per cm^2
M0 = 0.0;
for k = 1:Nz
    N_layer = N0_grid(k,:).';
    M0 = M0 + dz_cm * sum(N_layer .* v_vol);
end
fprintf('Initial mass in column   M0 = %.6e\n', M0);

%% --- RUN FULL SIMULATION WRAPPER ------------------------------------
sim   = CoagulationSimulation(cfg);
t_run = cfg.t_init:cfg.delta_t:cfg.t_final;

% tighter tolerances for better mass conservation
opts = odeset('RelTol',1e-8, 'AbsTol',1e-12, ...
              'NonNegative', 1:(Ns*Nz));

sim.run('tspan', t_run, 'v0', Y0, 'solver_options', opts);

t_all = sim.result.time;
Y_all = sim.result.concentrations;
Nt    = numel(t_all);

if Nt > 1
    dt_days = diff(t_all);
    dt_sec  = dt_days * 86400;
else
    dt_sec = [];
end

%% --- COLUMN MASS VS TIME --------------------------------------------
M_col = zeros(Nt,1);
for n = 1:Nt
    Yn   = Y_all(n,:).';
    Ymat = reshape(Yn, Ns, Nz).';   % [Nz x Ns]
    mass_n = 0.0;
    for k = 1:Nz
        N_layer = Ymat(k,:).';
        mass_n  = mass_n + dz_cm * sum(N_layer .* v_vol);
    end
    M_col(n) = mass_n;
end
Mf = M_col(end);

%% --- EXPORT AT 200 m (BOTTOM LAYER) ---------------------------------
r_cm = dg.getFractalRadii();
r_v  = dg.getConservedRadii();
ws   = SettlingVelocityService.velocity(r_cm, r_v, dg.setcon);  % [cm/s]
ws   = ws(:);

idx_start = (Nz-1)*Ns + 1;
idx_end   = Nz*Ns;

E_export = 0.0;
cumE     = zeros(Nt,1);

for n = 2:Nt
    layer_prev = Y_all(n-1, idx_start:idx_end).';
    layer_curr = Y_all(n,   idx_start:idx_end).';
    N_avg      = 0.5 * (layer_prev + layer_curr);

    flux_n = sum(N_avg .* v_vol .* ws);    % [biovol/cm^2/s]
    E_export   = E_export + flux_n * dt_sec(n-1);
    cumE(n)    = E_export;
end

fprintf('Final   mass in column   Mf = %.6e\n', Mf);
fprintf('Integrated export        E  = %.6e\n', E_export);
fprintf('Mf + E - M0                 = %.6e (should be ~0 if transport+coag conserve mass)\n', ...
    Mf + E_export - M0);
fprintf('Relative error              = %.3e\n', abs(Mf + E_export - M0)/M0);

%% --- PLOTS -----------------------------------------------------------
fig = figure('Color','w','Position',[100 100 1000 600]);

subplot(2,1,1);
plot(t_all, M_col, 'k', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Mass in Column');
title('1-D test: column mass vs time');
grid on;

subplot(2,1,2);
plot(t_all, cumE, 'b', 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Export (cum. mass)');
title('1-D test: export at 200 m vs time');
grid on;

outDir = fullfile(pwd, 'figures_tests');
if ~exist(outDir,'dir'), mkdir(outDir); end
print(fig, fullfile(outDir, 'Test1D_transport_only.png'), '-dpng','-r300');

fprintf('1-D transport(+coag) test complete.\n');
end