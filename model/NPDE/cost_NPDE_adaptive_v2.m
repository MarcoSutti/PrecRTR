function [ F ] = cost_NPDE_adaptive_v2( W0, W5_mat, pars )

% VERSIONE PRELIMINARE --- NOT THE MOST EFFICIENT IMPLEMENTATION!

% function [ F ] = cost_NPDE_adaptive_v2( W0, W5_mat, pars )
% Purpose: Computes the cost function value starting from the factorised
%          form of Wh.
% Created:     2022.10.19
% Last change: 2022.10.19

%   Oct 19, 2022:
%       Created by copying from cost_NPDE_f_from_struct.

global Gamma_h;
global Precomputed;

W0 = Dirichlet_on_struct( W0 );

G_1 = W0.U * W0.S;
H_1 = W0.V;
%--------------------------------------------------------------------------
n_h = size(W0.U,1);
h = 1/(n_h-1); 
area = h^2;
%--------------------------------------------------------------------------
lev = log2(n_h-1);
idx = lev - pars.lev_coarsest + 1;
%--------------------------------------------------------------------------
% Load from the global variable "Precomputed" the matrix L to calculate
% forward finite differences of the first derivatives
L = Precomputed(idx).L;
%--------------------------------------------------------------------------
G_Gamma = Gamma_h(idx).U * Gamma_h(idx).S;
H_Gamma = Gamma_h(idx).V;
%--------------------------------------------------------------------------
% Compute the factorized truncated elementwise product Wh.*Wh
% 04.12.2019: Keep the Wh_sq of the size k^2-by-k^2, i.e., do not truncate it!
% Wh_sq = getFactorizedHadamard( W0, W0 );
% 
% G_2 = Wh_sq.U * Wh_sq.S;
% H_2 = Wh_sq.V;
% 
% LG = L*G_1;
% LH = L*H_1;
% GtG = G_1'*G_1;
% HtH = H_1'*H_1;
% 
% first_term = 0.5 * trace( (LG'*LG) * HtH + GtG * (LH'*LH) );
% second_term = pars.lambda/2 * trace( GtG * HtH ) + pars.lambda/3 * trace( (G_1'*G_2) * (H_2'*H_1) ) ;
% third_term = - trace( (G_Gamma'*G_1) * (H_1'*H_Gamma) );
% 
% F = area * ( first_term + second_term + third_term );

W0_mat = W0.U * W0.S * W0.V';

W10 = W0_mat + W5_mat;

% Not efficient version, but fast to write :)
F = .5 * ( norm( L * W10, 'fro' )^2 + norm( W10 * L', 'fro' )^2 ) ...
    + .5 * pars.lambda * norm( W10, 'fro' )^2 + pars.lambda/3 * trace( W10' * W10.^2 ) ...
    - trace( H_Gamma * G_Gamma' * W10 );

F = area * F;

end