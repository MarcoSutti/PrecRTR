function [ EH ] = ehess_ACE( W, H, pars )

% function [ EH ] = ehess_ACE( W, H, pars )
% Purpose: Computes the directional derivative of the gradient along the
%          direction of H for the Allen-Cahn equation.
% Created:     2022.12.11
% Last change: 2023.04.11

%   Apr 11, 2023:
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

% Computes the Euclidean Hessian: 
sums_terms{1}  = pars.epsilon * pars.dt * ( A * H + H * A );
sums_terms{2} = (1-pars.dt) * H;
sums_terms{3} = 3*pars.dt * W.^2.*H;

EH = h_sq * sum(cat(3,sums_terms{:}),3);

end