function [gains, losses] = SectionalMassBalance(spec, p2)
%
% SectionalMassBalance calculates the sectional mass balance given a
% spectrum.
%
%

[n_times, n_sections] = size(spec);

% Calculate the gains & losses to each section. These come from coagulation
% and growth.

coag_gains  = zeros(n_times, n_sections);
coag_losses = zeros(n_times, n_sections);

for i_time = 1 : n_times
    
    vcon_r     = spec(i_time, :);
    vcon_shift = [0 vcon_r(1:n_sections-1)];
    
    term1 = vcon_r * p2.b2;
    term1 = vcon_r .* term1;

    term2 = vcon_r * p2.b1;
    term2 = term2 .* vcon_shift;
    
    coag_gains(i_time,:) = term1 + term2;
    
    term3 = vcon_r * (p2.b3 + p2.b4 + p2.b5);
    term3 = vcon_r .* term3;

    coag_losses(i_time, :) = term3;

end

sinking = diag(p2.sink_loss);
sinking = sinking';
sinking = sinking(ones(n_times,1),:);

sink_losses = sinking.*spec;

% Rearrange growth matrix (Note this may need to be changed)

g_gain    = diag(p2.growth,-1);
g_gain    = [p2.growth(1,1); g_gain];
g_gain    = g_gain';
g_gain    = g_gain(ones(i_time,1),:);
g_loss    = diag(p2.growth);
g_loss(1) = 0.0;
g_loss    = g_loss';
g_loss    = g_loss(ones(i_time,1),:);

growth_gain = g_gain.*spec;
growth_loss = g_loss.*spec;

gains.coag    = coag_gains;
gains.growth  = growth_gain;
losses.coag   = coag_losses;
losses.growth = growth_loss;
losses.settl  = sink_losses;
