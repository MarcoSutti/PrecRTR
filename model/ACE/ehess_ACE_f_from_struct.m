function [ EH ] = ehess_ACE_f_from_struct( W, H, pars )

% function [ EH ] = ehess_ACE_f_from_struct( W, H, pars )
% Purpose: Computes the directional derivative of the gradient along the
%          direction of H for the Allen-Cahn equation.
% Created:     2022.12.19
% Last change: 2023.01.11

%   Jan 11, 2023:
%       
%   Dec 19, 2022:
%       Created.

global Precomputed;

Nx = size(W.U,1);
Lx = 2*pi;
% Ly = 2*pi;
hx = Lx/Nx;
h_sq = hx^2;

%--------------------------------------------------------------------------
% Convert H to the ambient space format:
H_amb.U = [W.U*H.M + H.Up, W.U];
[ ~, k_h ] = size(H_amb.U);
H_amb.S = eye( k_h );
H_amb.V = [W.V, H.Vp];

%--------------------------------------------------------------------------
% W_sq = getFactorizedHadamard( W, W );
% W_sq_Hadamard_H = getFactorizedHadamard( W_sq, H_amb );

W_mat = W.U * W.S * W.V';
H_mat = H_amb.U * H_amb.S * H_amb.V';
W_sq_Hadamard_H_mat = (W_mat.^2).*H_mat;

% MS, 2023.02.03: This is very expensive! Besides the preconditioning 
%                 getXi, this is the most expensive line in the code for
%                 ACE.
[ Uk, Sk, Vk ] = getTruncatedSVD( W_sq_Hadamard_H_mat, pars.K );
W_sq_Hadamard_H.U = Uk;
W_sq_Hadamard_H.S = Sk;
W_sq_Hadamard_H.V = Vk;

%--------------------------------------------------------------------------
% Computes the Euclidean Hessian in factorized format:
% pref = -pars.epsilon * pars.dt * h_sq;

Ah = (1/h_sq) * Precomputed(1).Ah;

EH.U = [ ( pars.epsilon * pars.dt * Ah + (1-pars.dt) * speye(Nx) ) * H_amb.U, H_amb.U, W_sq_Hadamard_H.U ];
EH.S = h_sq * blkdiag( H_amb.S, pars.epsilon * pars.dt * H_amb.S, 3 * pars.dt * W_sq_Hadamard_H.S);
EH.V = [ H_amb.V, Ah * H_amb.V, W_sq_Hadamard_H.V ];


% % Check:
% W_mat = W.U * W.S * W.V';
% H_amb_mat = H_amb.U * H_amb.S * H_amb.V';
% EH_mat = EH.U * EH.S * EH.V';
% Hess_h = h_sq * (-pars.epsilon * pars.dt * D2 * H_amb_mat - pars.epsilon * pars.dt * H_amb_mat * D2 ...
%     + 3 * pars.dt * W_mat.^2 .* H_amb_mat + (1-pars.dt) * H_amb_mat );
% 
% % Hess_h - EH_mat

% norm( Hess_h - EH_mat, 'fro' )

end