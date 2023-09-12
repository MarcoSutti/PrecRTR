function [ EG ] = egrad_LYAP_adaptive_v2( W0, W5, pars )

% VERSIONE PRELIMINARE --- NOT THE MOST EFFICIENT IMPLEMENTATION!

% function [ EG ] = egrad_LYAP_adaptive_v2( W0, W5, pars )
% Purpose: Computes the Euclidean gradient in low-rank format.
% Created:     2022.10.06
% Last change: 2022.10.06

%   Oct 6, 2022:
%       Created by copying from egrad_LYAP_f_from_struct.

global Gamma_h;
global Precomputed;

W0 = Dirichlet_on_struct( W0 );

%--------------------------------------------------------------------------
[ n_h, ~ ] = size(W0.U);
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
EG.U = [ Ah * W0.U, W0.U, Gamma_h(idx).U, Ah * W5.U, W5.U ];
EG.V = [ W0.V, Ah * W0.V, Gamma_h(idx).V, W5.V, Ah * W5.V ];
EG.S = blkdiag( W0.S, W0.S, -area * Gamma_h(idx).S, W5.S, W5.S );

% Impose homogeneous Dirichlet BCs:
EG = Dirichlet_on_struct( EG );

end