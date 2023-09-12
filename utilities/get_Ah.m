function [ Ah ] = get_Ah( Nx )

% function [ Ah ] = get_Ah( Nx )
% Purpose: Compute the discretized minus Laplacian (without the prefactor
%          1/h^2) in sparse format, with periodic boundary conditions.
% Created:     2023.02.03
% Last change: 2023.02.03

%   Feb 3, 2023:
%       Created.

% The discretized minus Laplacian (without the prefactor 1/h^2) in sparse format:
Ah = spdiags( ones(Nx,1) * [-1 2 -1], -1:1, Nx, Nx );

% Impose periodic BCs on the discretized Laplacian:
Ah(1,end) = -1;
Ah(end,1) = -1;

% % Also possible construction using toeplitz(), but then the resulting
% matrix is not sparse.
% A2 = (1/hx^2) * toeplitz([-2, 1, zeros(1,Nx-3), 1]);

end