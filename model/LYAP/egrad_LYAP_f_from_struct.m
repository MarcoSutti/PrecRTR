function [ EG ] = egrad_LYAP_f_from_struct( Wh, pars )

% function [ EG ] = egrad_LYAP_f_from_struct( Wh, pars )
% Purpose: Computes the Euclidean gradient in low-rank format.
% Created:     29.05.2019
% Last change: 19.03.2020

%   March 19, 2020:
%       Added the global variable Precomputed to avoid recomputing Ah every
%       time the same level is visited.

global Gamma_h;
global Precomputed;

Wh = Dirichlet_on_struct( Wh );

%--------------------------------------------------------------------------
[ n_h, ~ ] = size(Wh.U);
h = 1/(n_h-1); 
area = h^2;
%--------------------------------------------------------------------------
lev = log2(n_h-1);
idx = lev - pars.lev_coarsest + 1;
%--------------------------------------------------------------------------
% Load the discretized minus Laplacian (without the prefactor 1/h^2) from
% the global variable "Precomputed"
Ah = Precomputed(idx).Ah;

% Computes the Euclidean gradient in factorized format:
EG.U = [ Ah * Wh.U, Wh.U, Gamma_h(idx).U ];
EG.V = [ Wh.V, Ah * Wh.V, Gamma_h(idx).V ];
EG.S = blkdiag( Wh.S, Wh.S, -area * Gamma_h(idx).S );

% Just to stay safe, impose homogeneous Dirichlet BCs:
EG = Dirichlet_on_struct( EG );

end