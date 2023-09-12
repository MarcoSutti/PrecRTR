function [ egrad_F ] = egrad_FKPP_CNLF( W, W0, W1, pars )

% function [ egrad_F ] = egrad_FKPP_CNLF( W, W0, W1, pars )
% Purpose: Computes the Euclidean gradient for the Fisher-KPP equation.
% Created:     2023.04.25
% Last change: 2023.04.25

%   Apr 25, 2023:
%       Created.

%--------------------------------------------------------------------------

star = pars.Mminus * W - pars.Mplus * W0 + 2*pars.dt*(W1.^2 - W1)*pars.Romega;

egrad_F = pars.Mminus'*star;

end