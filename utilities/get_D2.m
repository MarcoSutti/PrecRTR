function [ D2 ] = get_D2( L, N, h)

% function [ D2 ] = get_D2( L, N, h)
% Purpose: Calculate the second-order periodic spectral differentiation
% matrix of [Tre00, eq. (3.12)] rescaled to [−L, L].

% Reference: [Tre00] Trefethen, L. N. Spectral Methods in MATLAB. Society
%            for Industrial and Applied Mathematics, 2000.

% Created:     2022.12.09
% Last change: 2022.12.30

%   Dec 9, 2022:
%       Created.

column = [ -pi^2/(3*h^2)-1/6 -.5*(-1).^(1:N-1)./sin(h*(1:N-1)/2).^2 ];
D2 = (pi/L)^2 * toeplitz(column);  % 2022.12.30: Check this prefactor!!!!


% % 1) Get the second-order periodic spectral differentiation matrix:
% column = [ -pi^2/(3*hx^2)-1/6 -.5*(-1).^(1:Nx-1)./sin(hx*(1:Nx-1)/2).^2 ];
% D2 = (2*pi/Lx)^2 * toeplitz(column);   % 2022.12.30: avevo sbagliato il prefactor!!!

end