function [ Z ] = RK4( fun, Z0, h, Nt )

% function [ Z ] = RK4( fun, Z0, h, Nt )
% Purpose: Explicit Runge-Kutta method of 4th order.

% Reference: [] Wikipedia. Polycopié Analyse Numérique, p. ...

% fun: the right-hand side of the ODE.
% Z0:  the initial condition (a matrix).
% h:   the time step.
% Nt:  the total number of iterations.

% Output:
% Z:   a 3-way array containing the whole history of the RK4 time
%      integration.

% Created:     2022.12.29
% Last change: 2022.12.29

%   Dec 29, 2022:
%       Created by Marco Sutti.

Z(:,:,1) = Z0;

for n=1:Nt

    k1 = fun( Z(:,:,n) );
    k2 = fun( Z(:,:,n) + (h/2)*k1 );
    k3 = fun( Z(:,:,n) + (h/2)*k2 );
    k4 = fun( Z(:,:,n) + h*k3 );
    Z(:,:,n+1) = Z(:,:,n) + (h/6)*( k1 + 2*k2 + 2*k3 + k4 );
    
end
