function [ F ] = cost_FKPP_single_trajectory( w, w_tilde, pars )

% function [ F ] = cost_FKPP( w, w_tilde, pars )
% Purpose: Computes the cost function value for the Fisher-KPP equation.
% Created:     2023.04.11
% Last change: 2023.04.17

%   Apr 14, 2023:
%       New objective function.
%   Apr 11, 2023:
%       Created.

%--------------------------------------------------------------------------

% pref = .5 * pars.dt;
% 
% sums_terms(1) = pref * norm(pars.L * w, 'fro')^2;
% sums_terms(2) = .5 * sum(sum(w.^2));
% sums_terms(3) = -pref * sum(sum(w.^2 * pars.romega));
% sums_terms(4) = pars.dt/3 * sum(sum(w.^3 * pars.romega));
% sums_terms(5) = -sum(sum(w_tilde .* w_inner));
% 
% F = pars.hx * sum(sums_terms);

hr = pars.dt*pars.romega;
one_m_hr = 1 - hr;

% 2023.04.14: New functional!!!
F = .5 * norm( ( -pars.dt * pars.A + one_m_hr *speye(pars.Nx) ) * w + hr *w.^2 - w_tilde)^2;

end