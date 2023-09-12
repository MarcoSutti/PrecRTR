function [ EG ] = egrad_NPDE_adaptive_v2( W0, W5, pars )

% VERSIONE PRELIMINARE --- NOT THE MOST EFFICIENT IMPLEMENTATION!

% function [ EG ] = egrad_NPDE_adaptive_v2( W0, W5, pars )
% Purpose: Computes the Euclidean gradient in low-rank format.
% Created:     2022.10.19
% Last change: 2022.10.19

%   Oct 19, 2022:
%       Created by copying from cost_NPDE_f_from_struct.

global Gamma_h;
global Precomputed;

W0 = Dirichlet_on_struct( W0 );

%--------------------------------------------------------------------------
n_h = size(W0.U,1);
h = 1/(n_h-1); 
area = h^2;
%--------------------------------------------------------------------------
lev = log2(n_h-1);
idx = lev - pars.lev_coarsest + 1;
%--------------------------------------------------------------------------
% Load the discretized minus Laplacian (without the prefactor 1/h^2) from
% the global variable "Precomputed"
Ah = Precomputed(idx).Ah;
A_p_lambda_In = Ah + pars.lambda * area * speye(n_h);
%--------------------------------------------------------------------------
% 04.12.2019: Keep the Wh_squared of the size k^2-by-k^2, i.e., do not
% truncate it!
W0_sq = getFactorizedHadamard( W0, W0 );
W5_sq = getFactorizedHadamard( W5, W5 );
W_mix = getFactorizedHadamard( W5, W0 );
%--------------------------------------------------------------------------
% Computes the Euclidean gradient in factorized format:
EG.U = [ A_p_lambda_In * W5.U, W5.U, W5_sq.U, ...
         A_p_lambda_In * W0.U, W0.U, W0_sq.U, Gamma_h(idx).U, ...
         W_mix.U ];
EG.V = [ W5.V, Ah * W5.V, W5_sq.V, ...
         W0.V, Ah * W0.V, W0_sq.V, Gamma_h(idx).V, ... 
         W_mix.V ];
EG.S = blkdiag( W5.S, W5.S, pars.lambda * area * W5_sq.S, ...
                W0.S, W0.S, pars.lambda * area * W0_sq.S, -area * Gamma_h(idx).S, ...
                2 * pars.lambda * area * W_mix.S );

% Just to stay safe, impose homogeneous Dirichlet BCs:
EG = Dirichlet_on_struct( EG );

% % Computes the Euclidean gradient in factorized format:
% EG.U = [ Ah * W0.U, W0.U, Gamma_h(idx).U, Ah * W5.U, W5.U ];
% EG.V = [ W0.V, Ah * W0.V, Gamma_h(idx).V, W5.V, Ah * W5.V ];
% EG.S = blkdiag( W0.S, W0.S, -area * Gamma_h(idx).S, W5.S, W5.S );
% 
% % Impose homogeneous Dirichlet BCs:
% EG = Dirichlet_on_struct( EG );

end