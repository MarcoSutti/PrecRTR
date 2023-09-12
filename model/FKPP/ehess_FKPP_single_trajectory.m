function [ ehess_F ] = ehess_FKPP_single_trajectory( w, w_tilde, eta, pars )

% function [ ehess_F ] = ehess_FKPP_single_trajectory( w, w_tilde, H, pars )
% Purpose: Computes the directional derivative of the gradient along the
%          direction of H for the Fisher-KPP equation.
% Created:     2023.04.14
% Last change: 2023.04.17

%   Apr 15, 2023:
%       New Hessian due to new objective function.
%   Apr 14, 2023:
%       Created.

% %--------------------------------------------------------------------------
% % Computes the Euclidean Hessian: 
% sums_terms{1}  =  - pars.dt * pars.A * H;
% sums_terms{2} = H * (speye(pars.Nx) - pars.dt*pars.Romega);
% sums_terms{3} = 2 * pars.dt * w.*H * pars.Romega;
% 
% % Concatenate the cell array sums_terms along 3rd mode, then sum along the
% % same mode:
% EH = pars.hx * sum(cat(3,sums_terms{:}), 3);

hr = pars.dt*pars.romega;
one_m_hr = 1 - hr;
In = speye(pars.Nx);

% idem = one_m_hr * H - pars.dt * pars.A * H + 2 * hr * w.*H;
star = -pars.dt * pars.A * w + one_m_hr * w + hr *w.^2 - w_tilde;

primo_termine = 2*hr*(In.*eta)'*star;
first_ord_terms = - pars.dt * pars.A * eta + one_m_hr * eta + 2 * hr * w.*eta;
ehess_F = primo_termine + (- pars.dt * pars.A + one_m_hr * In + 2 * hr * (In.*w))'*first_ord_terms;

end