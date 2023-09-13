function [ EH ] = ehess_NPDE_f_from_struct( X, H, pars )

% function [ EH ] = ehess_NPDE_f_from_struct( X, H, pars )
% Purpose: Computes the directional derivative of the gradient along the
%          direction of H.
% Created:     10.03.2020
% Last change: 07.05.2020

%   May 7, 2020:
%       Changed the expression of the hessian so that the term \lambda W_{h}
%       is merged into AW_{h} and thereby reducing the rank of G_h.
%   March 19, 2020:
%       Added the global variable Precomputed to avoid recomputing Ah every
%       time the same level is visited.

global Precomputed;

%--------------------------------------------------------------------------
idx = pars.lev - pars.lev_coarsest + 1;
area = 2^(-2*pars.lev);
%--------------------------------------------------------------------------

X = Dirichlet_on_struct( X );

% Convert H to the ambient space format:
H_amb.U = [X.U*H.M + H.Up, X.U];
[ n_h, k_h ] = size(H_amb.U);
H_amb.S = speye( k_h );
H_amb.V = [X.V, H.Vp];

%--------------------------------------------------------------------------
% Load the discretized minus Laplacian (without the prefactor 1/h^2) from
% the global variable "Precomputed"
Ah = Precomputed(idx).Ah;
%--------------------------------------------------------------------------
X_Hadamard_H = getFactorizedHadamard( X, H_amb );
%--------------------------------------------------------------------------
% Computes the Euclidean Hessian in factorized format:
EH.U = [ ( Ah + pars.lambda * area * speye(n_h) ) * H_amb.U, H_amb.U, X_Hadamard_H.U ];
EH.V = [ H_amb.V, Ah * H_amb.V, X_Hadamard_H.V ];
EH.S = blkdiag( H_amb.S, H_amb.S, 2 * pars.lambda * area * X_Hadamard_H.S);

% Just to stay safe, impose homogeneous Dirichlet BCs:
EH = Dirichlet_on_struct( EH );

end