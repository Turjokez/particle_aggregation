% --- FILE: diag_biovolume_budget_run.m ---
function diag_biovolume_budget_run(sim, out)

t  = out.time(:);
od = out.output_data;
cfg = sim.config;

% authoritative inventory + export (already biovolume column units)
M_bv  = od.column_total_inventory_cm3cm2(:);   % cm^3/cm^2
EX_bv = od.bottom_total_flux_cm3cm2d(:);       % cm^3/cm^2/day

% unit flag
state_units = string(od.state_units);
fprintf('state_units = %s\n', state_units);

% RHS PP diagnostic
Y   = out.concentrations;
Ns  = cfg.n_sections;
Nz  = cfg.getNumLayers();
dzcm = cfg.dz * 100;

PP_bv = zeros(numel(t),1);

if state_units == "cm3"
    % state is already biovolume concentration: cm^3/cm^3
    for it = 1:numel(t)
        v = Y(it,:).';
        T = sim.rhs.decomposeTerms(t(it), v);
        dv_pp = reshape(T.dv_pp(:), [Ns Nz]);     % cm^3/cm^3/day
        PP_bv(it) = sum(dv_pp(:)) * dzcm;         % cm^3/cm^2/day
    end
else
    % state is number concentration: #/cm^3
    grd  = sim.grid;
    Vbin = grd.av_vol(:);                         % cm^3/particle
    conv = od.N_to_cm3;                           % model -> #/cm^3

    for it = 1:numel(t)
        v = Y(it,:).';
        T = sim.rhs.decomposeTerms(t(it), v);
        dv_pp = reshape(T.dv_pp(:), [Ns Nz]) * conv;  % #/cm^3/day
        PP_bv(it) = sum(sum(dv_pp .* Vbin, 1)) * dzcm; % cm^3/cm^2/day
    end
end

% integrated form budget
M_pred  = M_bv(1) + cumtrapz(t, PP_bv - EX_bv);
resid_M = M_bv - M_pred;

fprintf('max|resid| = %.3e (cm^3/cm^2)\n', max(abs(resid_M)));
fprintf('rel max|resid| = %.3e\n', max(abs(resid_M))/max(abs(M_bv)));

% plots
figure;
plot(t, M_bv, 'LineWidth',2); hold on; grid on;
plot(t, M_pred,'--','LineWidth',2);
legend('M_{bv} (model)','M_0 + \int(PP - EX)dt','Location','best');
xlabel('time (days)'); ylabel('cm^3/cm^2');
title('Biovolume inventory: model vs simple budget');

figure;
plot(t, resid_M, 'LineWidth',2); grid on;
xlabel('time (days)'); ylabel('cm^3/cm^2');
title('Simple biovolume budget residual');

end