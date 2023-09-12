function [ ehess_F ] = ehess_FKPP_CNLF_single_trajectory( w, eta, pars )

% function [ ehess_F ] = ehess_FKPP_CNLF_single_trajectory( w, eta, pars )
% Purpose: Computes the directional derivative of the gradient along the
%          direction of H for the Fisher-KPP equation.
% Created:     2023.04.24
% Last change: 2023.04.24

%   Apr 24, 2023:
%       Created.

%--------------------------------------------------------------------------

ehess_F = pars.MmtMm * eta;

end