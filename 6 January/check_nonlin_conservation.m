function check_nonlin_conservation(sim, t_target, z_target)

out = sim.result;
cfg = sim.config;
grd = sim.grid;

t = out.time(:);
[~, it] = min(abs(t - t_target));
tt = t(it);

Ns = cfg.n_sections;
Nz = cfg.getNumLayers();

% pick depth index closest to z_target
zc = grd.z_centers(:);
[~, kz] = min(abs(zc - z_target));

v = out.concentrations(it,:).';
terms = sim.rhs.decomposeTerms(tt, v);

% reshape to [Ns Nz]
dv_coag   = reshape(terms.dv_coag(:),   [Ns, Nz]);
dv_disagg = reshape(terms.dv_disagg(:), [Ns, Nz]);

% volume weights
if isprop(grd,'av_vol') && ~isempty(grd.av_vol)
    vol = grd.av_vol(:);
else
    r = grd.getConservedRadii();
    vol = (4/3)*pi*r(:).^3;
end

mc = sum(vol .* dv_coag(:,kz));
md = sum(vol .* dv_disagg(:,kz));
fprintf('\n=== Nonlinear conservation at t=%.4f d, z=%.1f m ===\n', tt, zc(kz));
fprintf('sum(vol * dv_coag)   = %+ .6e (should be ~0)\n', mc);
fprintf('sum(vol * dv_disagg) = %+ .6e (should be ~0 if disagg conserves same mass)\n', md);
fprintf('sum(vol * (coag+dis))= %+ .6e\n', mc+md);

end