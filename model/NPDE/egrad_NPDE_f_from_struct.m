function [ EG ] = egrad_NPDE_f_from_struct( Wh, pars )

% function [ EG ] = egrad_NPDE_f_from_struct( Wh, pars )
% Purpose: Computes the Euclidean gradient in low-rank format.
% Created:     23.07.2019
% Last change: 07.05.2020

%   May 7, 2020:
%       Changed the expression of the gradient so that the term \lambda W_{h}
%       is merged into AW_{h} and thereby reducing the rank of G_h.
%   March 19, 2020:
%       Added the global variable Precomputed to avoid recomputing Ah every
%       time the same level is visited.

global Gamma_h;
global Precomputed;

Wh = Dirichlet_on_struct( Wh );

%--------------------------------------------------------------------------
n_h = size(Wh.U,1);
h = 1/(n_h-1); 
area = h^2;
%--------------------------------------------------------------------------
lev = log2(n_h-1);
idx = lev - pars.lev_coarsest + 1;
%--------------------------------------------------------------------------
% Load the discretized minus Laplacian (without the prefactor 1/h^2) from
% the global variable "Precomputed"
Ah = Precomputed(idx).Ah;
%--------------------------------------------------------------------------
% 04.12.2019: Keep the Wh_squared of the size k^2-by-k^2, i.e., do not
% truncate it!
Wh_sq = getFactorizedHadamard( Wh, Wh );
%--------------------------------------------------------------------------
% Computes the Euclidean gradient in factorized format:
EG.U = [ ( Ah + pars.lambda * area * speye(n_h) ) * Wh.U, Wh.U, Wh_sq.U, Gamma_h(idx).U ];
EG.V = [ Wh.V, Ah * Wh.V, Wh_sq.V, Gamma_h(idx).V ];
EG.S = blkdiag( Wh.S, Wh.S, pars.lambda * area * Wh_sq.S, -area * Gamma_h(idx).S );

% Just to stay safe, impose homogeneous Dirichlet BCs:
EG = Dirichlet_on_struct( EG );

end