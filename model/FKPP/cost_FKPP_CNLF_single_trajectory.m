function [ F ] = cost_FKPP_CNLF_single_trajectory( w, w0, w1, pars )

% function [ F ] = cost_FKPP_CNLF_single_trajectory( w, w0, w1, pars )
% Purpose: Computes the cost function value for the Fisher-KPP equation.
% Created:     2023.04.24
% Last change: 2023.04.25

%   Apr 24, 2023:
%       Created.

%--------------------------------------------------------------------------

% star = pars.Mminus * w - pars.Mplus * w0 + 2 * pars.hr * (w1.^2 - w1);

% 2023.04.24: New functional!!!
% F = .5 * norm( star )^2;

Mmw = pars.Mminus * w;

F = .5 * (Mmw'*Mmw) - w0'*pars.Mplus'*Mmw + 2 * pars.hr * (w1.^2 - w1)' * Mmw;

end