function dvdt = CalcCoagDeriv(t, vcon, p2)
%
% CalcCoagDeriv calculates the derivates in the population balance
% equations
%
% Note, use a new, single structure to transfer all the parameters
%
% USAGE:
%
%
% HISTORY:
%  04-05-09: First cut - based heavily on GAJ's code
%
% Adrian Burd, University of Georgia, 2009

n_sections = length(vcon);

vcon(vcon<0) = eps;

vcon_r = vcon';      % vcon passed as column vector - make a row

vcon_shift = [0 vcon_r(1:n_sections-1)];

term1 = vcon_r * p2.b25;
term1 = vcon_r .* term1;

term2 = vcon_r * p2.b1;
term2 = term2 .* vcon_shift;

term3 = p2.linear * vcon;


%term4 = - diag(p2.disagg_minus).*vcon + [diag(p2.disagg_plus,-1);0].*[vcon(2:end);0];
%
%dvdt = (term1 + term2)' + term3 + term4;
dvdt = (term1 + term2)' + term3;

c3 = 0.2;
c4 = 1.45;

for isec = 2 : n_sections-1
    dvdt(isec) = dvdt(isec) - c3*c4^isec*(vcon(isec) - c4*vcon(isec+1));
end
