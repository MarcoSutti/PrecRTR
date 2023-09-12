function [ F ] = cost_LYAP_adaptive_v2( W0, W5_mat, pars )

% VERSIONE PRELIMINARE --- NOT THE MOST EFFICIENT IMPLEMENTATION!

% function [ F ] = cost_LYAP_adaptive_v2( W0, W5_mat, pars )
% Purpose: Computes the cost function value starting from the factorised
%          form of Wh.
% Created:     2022.10.06
% Last change: 2022.10.19

%   Oct 6, 2022:
%       Created by copying from cost_LYAP_f_from_struct.

global Gamma_h;
global Precomputed;

W0 = Dirichlet_on_struct( W0 );

% %--------------------------------------------------------------------------
% G = W0.U * W0.S;
% H = W0.V;
%--------------------------------------------------------------------------
n = size(W0.U,1);
h = 1/(n-1); 
area = h^2;
%--------------------------------------------------------------------------
lev = log2(n-1);
idx = lev - pars.lev_coarsest + 1;
%--------------------------------------------------------------------------
% Load from the global variable "Precomputed" the matrix L to calculate
% forward finite differences of the first derivative
L = Precomputed(idx).L;
%--------------------------------------------------------------------------
G_Gamma = Gamma_h(idx).U * Gamma_h(idx).S;
H_Gamma = Gamma_h(idx).V;
%--------------------------------------------------------------------------

% % See my notes of 24.06.2019
% sum_hadamard_prod = trace( (G_Gamma'*G) * (W0.V'*H_Gamma) );
% 
% % Factorized first derivatives
% LG = L*G;
% LH = L*H;
% 
% % forward difference
% F = 0.5 * area * trace( (LG'*LG) * (H'*H) + (LH'*LH) * (G'*G) ) ...
%     - area * sum_hadamard_prod;

W0_mat = W0.U * W0.S * W0.V';

W10 = W0_mat + W5_mat;

% Not efficient version, but fast to write :)
F = .5 * area * ( norm( L * W10, 'fro' )^2 + norm( W10 * L', 'fro' )^2 ) ...
    - area * trace( H_Gamma * G_Gamma' * W10 );

end