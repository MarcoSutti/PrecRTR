function [ EG ] = egrad_ACE_f_from_struct( W, W_tilde, pars )

% function [ EG ] = egrad_ACE_f_from_struct( W, W_tilde, pars )
% Purpose: Computes the Euclidean gradient in factorized low-rank format
%          for the Allen-Cahn equation.
% Created:     2022.12.19
% Last change: 2023.01.09

%   Dec 19, 2022:
%       Created.

global Precomputed;

Nx = size(W.U,1);
Lx = 2*pi;
% Ly = 2*pi;
hx = Lx/Nx;
h_sq = hx^2;
%--------------------------------------------------------------------------
% W_3rd Hadamard power
% W_3rd = getFactorizedHadamard( W, getFactorizedHadamard(W, W) );

W_mat = W.U * W.S * W.V';
 
[ Uk, Sk, Vk ] = svd( W_mat.^3 );
W_3rd.U = Uk;
W_3rd.S = Sk;
W_3rd.V = Vk;

%--------------------------------------------------------------------------
% Computes the Euclidean gradient in factorized format:
% pref = -pars.epsilon * pars.dt * h_sq;

Ah = (1/h_sq) * Precomputed(1).Ah;

EG.U = [ ( pars.epsilon * pars.dt * Ah + (1-pars.dt) * speye(Nx)) * W.U, W.U, W_3rd.U, W_tilde.U ];
EG.S = h_sq * blkdiag( W.S, pars.epsilon * pars.dt * W.S, pars.dt * W_3rd.S, - W_tilde.S );
EG.V = [ W.V, Ah * W.V, W_3rd.V, W_tilde.V ];

% % Load the discretized minus Laplacian (without the prefactor 1/h^2) from
% % the global variable "Precomputed"
% Ah = Precomputed(idx).Ah;
% 
% % Computes the Euclidean gradient in factorized format:
% EG.U = [ Ah * Wh.U, Wh.U, Gamma_h(idx).U ];
% EG.V = [ Wh.V, Ah * Wh.V, Gamma_h(idx).V ];
% EG.S = blkdiag( Wh.S, Wh.S, -area * Gamma_h(idx).S );

% % Check:
% W_mat = W.U * W.S * W.V';
% W_tilde_mat = W_tilde.U * W_tilde.S * W_tilde.V';
% EG_mat = EG.U * EG.S * EG.V';
% 
% first_term  = pref * D2 * W.U * W.S * W.V' + pref * W.U * W.S * W.V' * D2';
% 
% norm(EG_mat - first_term, "fro")
% pause

% second_term = (1-pars.dt) * h_sq * W.U * W.S * W.V';
% third_term = pars.dt * h_sq * W_3rd.U * W_3rd.S * W_3rd.V';
% fourth_term = - W_tilde.U * W_tilde.S * W_tilde.V';
% 
% first_term_mat  = pref * ( D2 * W_mat + W_mat * D2 );
% second_term_mat = (1-pars.dt) * h_sq * W_mat;
% third_term_mat = pars.dt * h_sq * W_mat.^3;
% fourth_term_mat = -W_tilde_mat;
% 
% % norm( first_term  - first_term_mat, 'fro')
% % norm( second_term - second_term_mat, 'fro')
% % norm( third_term  - third_term_mat, 'fro')
% norm( fourth_term - fourth_term_mat, 'fro' )
% pause

% % Computes the Euclidean gradient in factorized format:
% EG.U = [ Ah * Wh.U, Wh.U ];
% EG.V = [ Wh.V, Ah * Wh.V ];
% EG.S = blkdiag( Wh.S, Wh.S );

end