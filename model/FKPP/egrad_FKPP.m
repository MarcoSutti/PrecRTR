function [ egrad_F ] = egrad_FKPP( W, W_tilde, pars )

% function [ egrad_F ] = egrad_FKPP( W, W_tilde, pars )
% Purpose: Computes the Euclidean gradient for the Fisher-KPP equation.
% Created:     2023.04.11
% Last change: 2023.04.17

%   Apr 15, 2023:
%       New gradient due to new objective function.
%   Apr 11, 2023:
%       Created.

%--------------------------------------------------------------------------
% % Computes the Euclidean gradient: 
% sums_terms{1}  = -pars.dt * pars.A * W;
% sums_terms{2} = W * (speye(pars.Nx) - pars.dt*pars.Romega);
% sums_terms{3} = pars.dt * W.^2 * pars.Romega;
% sums_terms{4} = - W_tilde;
% 
% % Concatenate the cell array sums_terms along 3rd mode, then sum along the
% % same mode:
% egrad_F = pars.hx * sum(cat(3,sums_terms{:}), 3);

% W = W.U * W.S * W.V';

In = speye(pars.Nx);
hr = pars.dt * pars.Romega;
one_m_hr = In - hr;

star = -pars.dt * pars.A * W + W * one_m_hr + W.^2 * hr - W_tilde;

% diag_vecW = diag(W(:)); % spdiags(W(:),0,pars.Nx^2,pars.Nx^2);
diag_vecW = spdiags(W(:),0,pars.Nx^2,pars.Nx^2);
star_Romega = star * pars.Romega;
third_term_vec = diag_vecW * star_Romega(:);

third_term = 2 * pars.dt*reshape(third_term_vec,[pars.Nx,pars.Nx]);

egrad_F = -pars.dt * pars.A'*star + star*one_m_hr + third_term;

end