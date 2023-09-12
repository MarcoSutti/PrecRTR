function [ F ] = cost_FKPP_CNLF( W, W0, W1, pars )

% function [ F ] = cost_FKPP_CNLF( w, w0, w1, pars )
% Purpose: Computes the cost function value for the Fisher-KPP equation.
% Created:     2023.04.25
% Last change: 2023.04.25

%   Apr 25, 2023:
%       Created.

%--------------------------------------------------------------------------

star = pars.Mminus * W - pars.Mplus * W0 + 2 * pars.dt * (W1.^2 - W1) * pars.Romega;

F = .5 * norm( star, 'fro')^2;

% F = .5 * (Mmw'*Mmw) - W0'*pars.Mplus'*Mmw + 2 * pars.hr * (W1.^2 - W1)' * Mmw;

end