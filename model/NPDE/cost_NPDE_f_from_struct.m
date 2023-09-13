function [ F ] = cost_NPDE_f_from_struct( Wh, pars )

% function [ F ] = cost_NPDE_f_from_struct( Wh, pars )
% Purpose: Computes the cost function value starting from the factorised
%          form of Wh.
% Created:     23.07.2019
% Last change: 31.03.2020

%   March 19, 2020:
%       Added the global variable Precomputed to avoid recomputing Ah every
%       time the same level is visited.

global Gamma_h;
global Precomputed;

Wh = Dirichlet_on_struct( Wh );

G_1 = Wh.U * Wh.S;
H_1 = Wh.V;
%--------------------------------------------------------------------------
n_h = size(Wh.U,1);
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
Wh_sq = getFactorizedHadamard( Wh, Wh );

G_2 = Wh_sq.U * Wh_sq.S;
H_2 = Wh_sq.V;

LG = L*G_1;
LH = L*H_1;
GtG = G_1'*G_1;
HtH = H_1'*H_1;

first_term = 0.5 * trace( (LG'*LG) * HtH + GtG * (LH'*LH) );
second_term = pars.lambda/2 * trace( GtG * HtH ) + pars.lambda/3 * trace( (G_1'*G_2) * (H_2'*H_1) ) ;
third_term = - trace( (G_Gamma'*G_1) * (H_1'*H_Gamma) );

F = area * ( first_term + second_term + third_term );

end