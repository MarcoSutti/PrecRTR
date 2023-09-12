function [ egrad_F ] = egrad_FKPP_CNLF_single_trajectory( w, w0, w1, pars )

% function [ egrad_F ] = egrad_FKPP_CNLF_single_trajectory( w, w0, w1, pars )
% Purpose: Computes the Euclidean gradient for the Fisher-KPP equation.
% Created:     2023.04.24
% Last change: 2023.04.25

%   Apr 24, 2023:
%       Created.

%--------------------------------------------------------------------------

star = pars.Mminus * w - pars.Mplus * w0 + 2*pars.hr*(w1.^2 - w1);

egrad_F = pars.Mminus'*star;

end