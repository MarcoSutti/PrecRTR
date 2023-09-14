function [ Ah ] = get_Ah_Neumann_BCs( Nx, hx )

% function [ Ah ] = get_Ah_Neumann_BCs( Nx, hx )
% Purpose: Compute the discretized Laplacian in sparse format, with Neumann
%          boundary condition encoded.
% Created:     2023.04.19
% Last change: 2023.04.27

%   Apr 27, 2023:
%       Changed the sign, and added the prefactor (1/hx^2).
%   Apr 19, 2023:
%       Created.

% The discretized Laplacian (without the prefactor 1/h^2) in sparse format:
Ah = spdiags( ones(Nx,1) * [1 -2 1], -1:1, Nx, Nx );

% Impose Neumann BCs on the discretized Laplacian:
Ah(end,end-1) = 2;
Ah(1,2) = 2;

% Add the prefactor 1/h^2:
Ah = (1/hx^2) * Ah;

end