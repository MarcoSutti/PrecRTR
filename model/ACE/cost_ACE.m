function [ F ] = cost_ACE( W, W_tilde, pars )

% function [ F ] = cost_ACE( W, W_tilde, pars )
% Purpose: Computes the cost function value for the Allen-Cahn equation.
% Created:     2023.01.23
% Last change: 2023.04.11

%   Jan 23, 2023:
%       Created.

%--------------------------------------------------------------------------
Nx = size(W,1);
Lx = 2*pi;
% Ly = 2*pi;
hx = Lx/Nx;
h_sq = hx^2;
%--------------------------------------------------------------------------

% Create matrix L in sparse format
% The matrix L is used to calculate forward finite differences of the
% first derivatives
B = [ -ones(Nx,1), [0; ones(Nx-1,1)] ];
d = [0,1];
L = (1/hx)*spdiags( B, d, Nx, Nx );

% Impose periodic BCs on L:
L(end,1) = L(1,2);
%--------------------------------------------------------------------------

pref = .5 * pars.epsilon * pars.dt;

sums_terms(1) = pref * ( norm(L * W, 'fro')^2 + norm(W * L', 'fro')^2 );

sums_terms(2) = .5 * (1-pars.dt) * sum(sum(W.^2));

sums_terms(3) = .25 * pars.dt * sum(sum(W.^4));

sums_terms(4) = -sum(sum(W_tilde.*W));
% sums_terms(4) = - trace( W_tilde' * W );
 
F = h_sq * sum(sums_terms);

end