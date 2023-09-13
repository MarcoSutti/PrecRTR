function [ EH ] = ehess_LYAP_f_from_struct( X, H, pars )

% function [ EH ] = ehess_LYAP_f_from_struct( X, H, pars )
% Purpose: Computes the Euclidean Hessian in factorized format.
% Created:     10.03.2020
% Last change: 19.03.2020

%   March 19, 2020:
%       Added the global variable Precomputed to avoid recomputing Ah every
%       time the same level is visited.

global Precomputed;

%--------------------------------------------------------------------------
idx = pars.lev - pars.lev_coarsest + 1;
%--------------------------------------------------------------------------

X = Dirichlet_on_struct( X );

% Convert H to the ambient space format:
H_amb.U = [X.U*H.M + H.Up, X.U];
k_h = size(H_amb.U,2);
H_amb.S = speye( k_h );
H_amb.V = [X.V, H.Vp];

%--------------------------------------------------------------------------
% Load the discretized minus Laplacian (without the prefactor 1/h^2) from
% the global variable "Precomputed"
Ah = Precomputed(idx).Ah;

% Computes the Euclidean Hessian in factorized format:
EH.U = [ Ah * H_amb.U, H_amb.U ];
EH.V = [ H_amb.V, Ah * H_amb.V ];
EH.S = blkdiag( H_amb.S, H_amb.S );

% Just to stay safe, impose homogeneous Dirichlet BCs:
EH = Dirichlet_on_struct( EH );

end