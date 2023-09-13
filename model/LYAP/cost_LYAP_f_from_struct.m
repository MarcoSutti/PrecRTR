function [ F ] = cost_LYAP_f_from_struct( Wh, pars )

% function [ F ] = cost_LYAP_f_from_struct( Wh, pars )
% Purpose: Computes the cost function value starting from the factorised
%          form of Wh.
% Created:     17.04.2019
% Last change: 19.03.2020

%   March 19, 2020:
%       Added the global variable Precomputed to avoid recomputing Ah every
%       time the same level is visited.

global Gamma_h;
global Precomputed;

Wh = Dirichlet_on_struct( Wh );

%--------------------------------------------------------------------------
G = Wh.U * Wh.S;
H = Wh.V;
%--------------------------------------------------------------------------
n_h = size(Wh.U,1);
h = 1/(n_h-1); 
area = h^2;
%--------------------------------------------------------------------------
lev = log2(n_h-1);
idx = lev - pars.lev_coarsest + 1;
%--------------------------------------------------------------------------
% Load from the global variable "Precomputed" the matrix L to calculate
% forward finite differences of the first derivative
L = Precomputed(idx).L;
%--------------------------------------------------------------------------
G_Gamma = Gamma_h(idx).U * Gamma_h(idx).S;
H_Gamma = Gamma_h(idx).V;
%--------------------------------------------------------------------------

% See my notes of 24.06.2019
sum_hadamard_prod = trace( (G_Gamma'*G) * (Wh.V'*H_Gamma) );

% Factorized first derivatives
LG = L*G;
LH = L*H;

% forward difference
F = 0.5 * area * trace( (LG'*LG) * (H'*H) + (LH'*LH) * (G'*G) ) ...
    - area * sum_hadamard_prod;

end