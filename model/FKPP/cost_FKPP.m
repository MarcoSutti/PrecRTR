function [ F ] = cost_FKPP( W, W_tilde, pars )

% function [ F ] = cost_FKPP( W, W_tilde, pars )
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
% sums_terms(1) = pref * norm(pars.L * W, 'fro')^2;
% sums_terms(2) = .5 * sum(sum(W.^2));
% sums_terms(3) = -pref * sum(sum(W.^2 * pars.Romega));
% sums_terms(4) = pars.dt/3 * sum(sum(W.^3 * pars.Romega));
% sums_terms(5) = -sum(sum(W_tilde .* W));
% 
% F = pars.hx * sum(sums_terms);

hr = pars.dt * pars.Romega;
one_m_hr = speye(pars.Nx) - hr;

% 2023.04.14: New functional!!!
F = .5 * norm( -pars.dt * pars.A * W + W * one_m_hr + W.^2 * hr - W_tilde, 'fro')^2;

end