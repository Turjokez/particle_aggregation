function outflag = CoagOutput(p, p2, t_out, spec)
%
% This function does some simple graphical output
%

outflag = 0;

%% Diagnostics
%  Run some diagnostics on the model output

[sec_gains, sec_losses]     = SectionalMassBalance(spec, p2);
[total_gains, total_losses] = TotalMassBalance(spec, p2);

%% Book-keeping
%  Here we calculate some things for the graphical output

total_spec = sum(spec');

axis_limits = [0, p.n_sections, min(total_spec), max(total_spec)];

n_times = length(t_out);

%% Calculate additional spectra
%  Calculate the number spectrum, mass spectrum and flux with respect to
%  particle diameters.

nspec_v    = zeros(n_times, p.n_sections);
massspec_v = nspec_v;
fluxsect   = nspec_v;
fluxspec   = nspec_v;

r_i = p.amfrac *p.av_vol.^p.bmfrac;
r_v = (0.75/pi*p.av_vol).^(1.0/3.0);

set_vel = SettlingVelocity(r_i, r_v, p.setcon);
set_vel = set_vel/100*p.day_to_sec;

diam_i = 2.0*p.r_to_rg*r_i;
diam_v = 2.0*r_v;

diam_i = diam_i';
diam_v = diam_v';

diam_i_mat = diam_i(ones(n_times, 1), :);
diam_v_mat = diam_v(ones(n_times, 1), :);

for jindx = 1 : n_times

    yout = spec(jindx,:);
    nspec_v(jindx,:)   = yout./(1.5*p.v_lower')./p.dwidth';
    masspec_v(jindx,:) = yout./p.dwidth';
    fluxsect(jindx,:)  = yout.*set_vel'*1e6;
    fluxspec(jindx,:)  = masspec_v(jindx,:).*set_vel'*1e6;

end
total_flux = sum(fluxsect,2);
total_mass = sum(spec, 2);

diaratio = (p.fr_dim/3)*diam_v_mat./diam_i_mat;
nspec_i = nspec_v.*diaratio;
masspec_i = masspec_v.*diaratio;
fluxspec_i = fluxspec.*diaratio;


%% Simple 2D plots

fig_h1 = figure(1);

sp_h_221 = subplot(2,2,1);
loglog(diam_i, nspec_i(1,:), 'b', diam_i, nspec_i(end,:), 'r');
xlabel('Particle diameter [cm]')
ylabel('Number spectrum [# cm^{-4}]')
axis tight

sp_h_223 = subplot(2,2,3);
plot(diam_i_mat(:,2:end)', fluxspec_i(:,2:end)')
xlabel('Particle image diameter [cm]')
ylabel('Volume flux spectra [cm^2 m^{-2} d^{-1}]')
axis tight

sp_h_247 = subplot(2,4,7);
plot(t_out, fluxsect, t_out, total_flux, '*--')
xlabel('Time [d]')
ylabel('Sectional Flux [cm^3 m^{-2} d^{-1} sect^{-1}]')

sp_h_248 = subplot(2,4,8);
plot(t_out, total_flux./total_mass/1e6)
xlabel('Time [d]')
ylabel('Average v [m d^{-1}]')

sp_h_243 = subplot(2,4,3);
semilogy(t_out, spec, t_out, total_mass, '*--')
xlabel('Time [d]')
ylabel('Sectional concentration [vol/vol/sect]')

fig_h2 = figure(2);
plot(t_out, total_gains.growth./(total_losses.sett + total_losses.coag))
xlabel('Time [d]')
ylabel('Gains/Losses')
title('Total System Mass Balance')

fig_h3 = figure(3);
plot(t_out, total_losses.coag./total_losses.sett)
set(gca, 'FontName', 'Helvetica', 'FontSize', 14)
xlabel('Time [d]', 'FontName', 'Helvetica', 'FontSize', 14)
ylabel('(Coag Losses)/(Settling Losses)', 'FontName', 'Helvetica', 'FontSize', 14)

t_mat = t_out(:,ones(1,p.n_sections));
v_mat = p.v_lower';
v_mat = v_mat(ones(length(t_out),1),:);

% %keyboard
%
% fig_h4 = figure(4);
% surf(log10(v_mat), t_mat, log10(sec_gains.coag))
% xlabel('Volume')
% ylabel('Time [d]')
% zlabel('Coagulation Gains')
%
% fig_h5 = figure(5);
% surf(log10(v_mat), t_mat, log10(sec_losses.coag))
% xlabel('Volume')
% ylabel('Time [d]')
% zlabel('Coagulation Losses')

sectional_gains = sec_gains.coag + sec_gains.growth;
sectional_loss  = sec_losses.coag + sec_losses.settl + sec_losses.growth;

% Debug: Print legacy data ranges for comparison
fprintf('LEGACY DEBUG: sectional_gains range: [%.2e, %.2e]\n', min(sectional_gains(:)), max(sectional_gains(:)));
fprintf('LEGACY DEBUG: sectional_loss range: [%.2e, %.2e]\n', min(sectional_loss(:)), max(sectional_loss(:)));
fprintf('LEGACY DEBUG: ratio range: [%.2e, %.2e]\n', min(sectional_gains(:)./sectional_loss(:)), max(sectional_gains(:)./sectional_loss(:)));
fprintf('LEGACY DEBUG: v_mat range: [%.2e, %.2e]\n', min(v_mat(:)), max(v_mat(:)));

fig_h6 = figure(6);
surf(log10(v_mat), t_mat, sectional_gains./sectional_loss)
xlabel('Volume')
ylabel('Time [d]')
zlabel('Gains/Losses')



%% Save diagnostics for comparison with OOP
try
    diag_legacy = struct();
    diag_legacy.t = t_out;
    diag_legacy.Y = spec;
    diag_legacy.r_i = r_i;
    diag_legacy.r_v = r_v;
    diag_legacy.diam_i = diam_i;
    diag_legacy.diam_v = diam_v;
    diag_legacy.set_vel = set_vel;
    diag_legacy.v_lower = p.v_lower;
    diag_legacy.dwidth = p.dwidth;
    diag_legacy.nspec_v = nspec_v;
    diag_legacy.masspec_v = masspec_v;
    diag_legacy.fluxspec = fluxspec;
    diag_legacy.diaratio = diaratio;
    diag_legacy.fluxspec_i = fluxspec_i;
    diag_legacy.total_losses = total_losses;
    diag_legacy.total_gains = total_gains;
    save('plot_diag_legacy.mat','diag_legacy');
catch
end