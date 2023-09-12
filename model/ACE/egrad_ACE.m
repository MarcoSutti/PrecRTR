function [ egrad_F ] = egrad_ACE( W, W_tilde, pars )

% function [ egrad_F ] = egrad_ACE( W, W_tilde, pars )
% Purpose: Computes the Euclidean gradient for the Allen-Cahn equation.
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
% The discretized minus Laplacian (without the prefactor 1/h^2) in sparse format:
A = spdiags( ones(Nx,1) * [-1 2 -1], -1:1, Nx, Nx );

% Impose periodic BCs on the discretized Laplacian:
A(1,end) = -1;
A(end,1) = -1;

A = (1/h_sq) * A;
%--------------------------------------------------------------------------

% Computes the Euclidean gradient: 
sums_terms{1}  = pars.epsilon * pars.dt * ( A * W + W * A );
sums_terms{2} = (1-pars.dt) * W;
sums_terms{3} = pars.dt * W.^3;
sums_terms{4} = - W_tilde;

egrad_F = h_sq * sum(cat(3,sums_terms{:}),3);

end