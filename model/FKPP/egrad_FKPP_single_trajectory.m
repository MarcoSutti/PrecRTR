function [ egrad_F ] = egrad_FKPP_single_trajectory( w, w_tilde, pars )

% function [ egrad_F ] = egrad_FKPP_single_trajectory( W, W_tilde, pars )
% Purpose: Computes the Euclidean gradient for the Fisher-KPP equation.
% Created:     2023.04.14
% Last change: 2023.04.17

%   Apr 15, 2023:
%       New gradient due to new objective function.
%   Apr 14, 2023:
%       Created.

%--------------------------------------------------------------------------
% % Computes the Euclidean gradient: 
% sums_terms{1}  = -pars.dt * pars.A * w;
% sums_terms{2} = (1 - pars.dt*pars.romega) * w;
% sums_terms{3} = pars.dt * w.^2 * pars.romega;
% sums_terms{4} = - w_tilde;
% 
% % Concatenate the cell array sums_terms along 3rd mode, then sum along the
% % same mode:
% egrad_F = pars.hx * sum(cat(3,sums_terms{:}), 3);

hr = pars.dt*pars.romega;
one_m_hr = 1 - hr;

calcolato_una_volta = -pars.dt * pars.A * w + one_m_hr * w + hr *w.^2 - w_tilde;

egrad_F = (one_m_hr*speye(pars.Nx) - pars.dt*pars.A')*calcolato_una_volta + 2*hr *calcolato_una_volta.*w;

end