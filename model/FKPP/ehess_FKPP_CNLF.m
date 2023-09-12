function [ ehess_F ] = ehess_FKPP_CNLF( W, H, pars )

% function [ ehess_F ] = ehess_FKPP_CNLF( W, H, pars )
% Purpose: Computes the directional derivative of the gradient along the
%          direction of H for the Fisher-KPP equation.
% Created:     2023.04.25
% Last change: 2023.04.25

%   Apr 25, 2023:
%       Created.

%--------------------------------------------------------------------------

ehess_F = pars.MmtMm * H;

end