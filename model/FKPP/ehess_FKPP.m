function [ EH ] = ehess_FKPP( H, pars )

% function [ EH ] = ehess_FKPP( H, pars )
% Purpose: Computes the directional derivative of the gradient along the
%          direction of H for the Fisher-KPP equation.
% Created:     2023.04.11
% Last change: 2023.05.08

%   May 8, 2023:
%       Changing the cost function changed everything.
%   Apr 15, 2023:
%       New Hessian due to new objective function.
%   Apr 11, 2023:
%       Created.

%--------------------------------------------------------------------------
% % Computes the Euclidean Hessian: 
% sums_terms{1}  =  - pars.dt * pars.A * H;
% sums_terms{2} = H * (speye(pars.Nx) - pars.dt*pars.Romega);
% sums_terms{3} = 2 * pars.dt * W.*H * pars.Romega;
% 
% % Concatenate the cell array sums_terms along 3rd mode, then sum along the
% % same mode:
% EH = pars.hx * sum(cat(3,sums_terms{:}), 3);

% In = speye(pars.Nx);
% hr = pars.dt * pars.Romega;
% one_m_hr = In - hr;
% 
% star = -pars.dt * pars.A * W + W * one_m_hr + W.^2 * hr - W_tilde;
% square = -pars.dt * pars.A * H + H * one_m_hr + 2*(W.*H)*hr;
% 
% diag_vecH = diag(H(:)); % spdiags(W(:),0,pars.Nx^2,pars.Nx^2);
% star_Romega = star * pars.Romega;
% first_term_vec = diag_vecH * star_Romega(:);
% first_term = 2 * pars.dt *reshape(first_term_vec,[pars.Nx,pars.Nx]);
% 
% diag_vecW = diag(W(:)); % spdiags(W(:),0,pars.Nx^2,pars.Nx^2);
% square_Romega = square * pars.Romega;
% fourth_term_vec = diag_vecW * square_Romega(:);
% fourth_term = 2 * pars.dt *reshape(fourth_term_vec,[pars.Nx,pars.Nx]);
% 
% EH = first_term - pars.dt*pars.A'*square + square*one_m_hr + fourth_term;

EH = pars.MmtMm * H;

end