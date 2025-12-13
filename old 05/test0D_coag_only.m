function test0D_coag_only()
% TEST0D_COAG_ONLY
% Pure 0-D coagulation test (no growth, no sinking, no breakup, no losses).
% Checks that total biovolume is conserved in time.

clc;
fprintf('=============================================\n');
fprintf('   TEST 0-D: COAGULATION ONLY\n');
fprintf('=============================================\n');

%% --- CONFIG (SLAB ONLY) ----------------------------------------------
cfg = SimulationConfig();
cfg.use_column   = false;        % 0-D slab
cfg.n_sections   = 35;

cfg.t_init       = 0.0;
cfg.t_final      = 2.0;          % [days]
cfg.delta_t      = 0.1;

% physics for a *gentle* coagulation
cfg.epsilon_profile = 'constant';
cfg.epsilon         = 1e-9;
cfg.alpha_base      = 0.01;

% NO BIOLOGY / LOSSES
cfg.use_NPP                  = false;
cfg.growth                   = 0.0;
cfg.attenuation_rate         = 0.0;
cfg.attenuation_depth_factor = 0.0;
cfg.loss_rate                = 0.0;

% NO NONLINEAR BREAKUP
cfg.disagg_use_nonlinear     = false;

%% --- GRID + BETA MATRICES -------------------------------------------
dg    = cfg.derive();
Ns    = cfg.n_sections;
v_vol = dg.av_vol(:);      % [cm^3]

% Brownian + shear + diff. sedimentation kernels
assembler = BetaAssembler(cfg, dg);
b_brown   = assembler.computeFor('KernelBrown');
b_shear   = assembler.computeFor('KernelCurSh');
b_ds      = assembler.computeFor('KernelCurDS');
betas     = assembler.combineAndScale(b_brown, b_shear, b_ds);

%% --- PURE COAGULATION RHS -------------------------------------------
linear_zero = zeros(Ns);          % no growth / sinking
dis0        = zeros(Ns);          % no linear disaggregation

rhs = CoagulationRHS( ...
    betas, ...
    linear_zero, ...  % linear
    dis0, ...         % disaggMinus
    dis0, ...         % disaggPlus
    cfg);

%% --- INITIAL CONDITION ----------------------------------------------
% simple: 1e3 particles in each section
N0 = 1e3 * ones(Ns,1);

M0 = sum(N0 .* v_vol);
fprintf('Initial mass M0  = %.6e\n', M0);

%% --- INTEGRATE USING ODE15S DIRECTLY --------------------------------
t_grid = cfg.t_init : cfg.delta_t : cfg.t_final;

% use tight tolerances + non-negative constraint
opts = odeset('RelTol',1e-8, 'AbsTol',1e-12, ...
              'NonNegative',1:Ns);

% rhs function for ode15s
rhsFun = @(t,y) rhs.evaluate(t,y);

fprintf('Integrating 0-D coagulation (%d steps)...\n', numel(t_grid));
[t_out, Y_all] = ode15s(rhsFun, t_grid, N0, opts);

%% --- MASS DIAGNOSTIC -----------------------------------------------
M_t = Y_all * v_vol;      % [Nt x 1]
Mf  = M_t(end);

fprintf('Final   mass Mf  = %.6e\n', Mf);
fprintf('Mf - M0          = %.6e (should be ~0 if coag is conserving)\n', Mf - M0);
fprintf('Relative error   = %.3e\n', abs(Mf-M0)/M0);

%% --- PLOT -----------------------------------------------------------
fig = figure('Color','w','Position',[100 100 900 550]);
plot(t_out, M_t, 'LineWidth', 2);
xlabel('Time (days)');
ylabel('Total mass (biovolume)');
title('0-D Coagulation Only: Mass vs Time');
grid on;

outDir = fullfile(pwd, 'figures_tests');
if ~exist(outDir,'dir'), mkdir(outDir); end
print(fig, fullfile(outDir, 'Test0D_coag_only.png'), '-dpng','-r300');

fprintf('0-D coagulation test complete.\n');
end