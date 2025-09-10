function dfdy = CalcCoagJac(t, vcon, p2)
%
% CalcCoagJac calculates the Jacobian for the coagulation equations. Each
% row corresponds to a different function and each column to the derivative
% of that function with respect to a different variable. So, if
%
%      f_i(y) = dQ_i/dt
%
%    J_{ij} = (df_1/dy_1   df_1/dy_2  ... df_1/dy_n)
%             (df_2/dy_1   df_2/dy_2  ... df_2/dy_n)
%             (   :            :              :    )
%             (   :            :              :    )
%             (df_m/dy_1   df_m/dy_2  ... df_m/dy_n)
%
%
% USAGE
%
% HISTORY
%   04-05-09: First cut
%
% Adrian Burd, University of Georgia
%
n_sections = length(vcon);

vcon_r = vcon';

vcon_mat = vcon_r(ones(1, n_sections),:);
vcon_shift = [zeros(n_sections,1) vcon_mat(:,1:end-1)];

% First calculate the df_i/dy_i terms

term1 = vcon_r*p2.b25;
term1 = diag(term1);

term1 = term1 + diag(vcon_r).*p2.b25;

% Calculate the df_i/dy_{i-1} terms

term2a = vcon_r*p2.b1;
term2a = diag(term2a(2:end),-1);

term2b = diag(p2.b1, 1);
term2b = term2b' .* vcon_r(1:end-1);
term2b = diag(term2b, -1);

% Change the following line.
%term2c = diag(vcon_r(1:end-1),-1) .* p2.b25;

term2c = diag(vcon_r(2:end),-1) .* p2.b25';

term2 = term2a + term2b + term2c;

% Calculate the df_i/dy_j terms (j ~= i, i-1)

% term3a = triu(p2.b1, 2) .* vcon_shift;
% term3a = term3a';
%
% term3b = (triu(p2.b25, 2) + tril(p2.b25, -1)) .* vcon_mat;
%
% term3 = term3b';    %%% CHECK THIS, WHY NO TERM3A?

term3a = p2.b1 .* vcon_shift;

term3b = p2.b25 .* vcon_mat;

term3 = (term3a + term3b)';

term3 = triu(term3,2) + tril(term3,-1);

% Now the linear terms

lin_term = p2.linear;

dfdy = term1 + term2 + term3 + lin_term - p2.disagg_minus + p2.disagg_plus;

% % Test using loops
%
% dfdy2a = zeros(n_sections, n_sections);
% dfdy2b = zeros(n_sections, n_sections);
% dfdy2c = zeros(n_sections, n_sections);
% dfdy2ca = zeros(n_sections, n_sections);
% dfdy2cb = zeros(n_sections, n_sections);
%
% sum1 = 0.0;
% sum2 = 0.0;
%
% for isec = 1 : n_sections
%     sum1 = 0.0;
%     for ksec = 1 : n_sections
%         sum1 = sum1 + p2.b25(ksec, isec)*vcon(ksec);
%     end
%     dfdy2a(isec,isec) = sum1 + p2.b25(isec,isec)*vcon(isec);
% end
%
% for isec = 2 : n_sections
%     sum2 = 0.0;
%     if(isec > 1)
%         for ksec = 1 : isec -1
%             sum2 = sum2 + p2.b1(ksec,isec)*vcon(ksec);
%         end
%     end
%     dfdy2b(isec,isec-1) = sum2 + p2.b1(isec-1,isec)*vcon(isec-1) + p2.b25(isec-1,isec)*vcon(isec);
% end
%
% for isec = 1 : n_sections
%    for jsec = 1 : n_sections
%        if (jsec ~= isec || jsec ~= isec-1)
%           dfdy2ca(isec,jsec) = p2.b25(jsec,isec)*vcon(isec);
%           if (isec > 1)
%              dfdy2ca(isec,jsec) = p2.b1(jsec,isec)*vcon(isec-1);
%           end
%           [isec jsec dfdy2ca(isec,jsec) dfdy2ca(isec,jsec)]
%
%           dfdy2c(isec,jsec) = dfdy2ca(isec,jsec) + dfdy2cb(isec,jsec);
%        end
%    end
% end
%
% keyboard