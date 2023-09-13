function [ Gamma ] = RHS_NPDE_struct( h )

% function [ Gamma ] = RHS_NPDE_struct( h )
% Purpose: Computes the right-hand side for the NPDE problem.
% Created:     30.07.2019
% Last change: 30.03.2020

%   March 30, 2020:
%       Changed the name of this function to "RHS_NPDE_struct".

x = (0:h:1)';
y = (0:h:1)';

% A very simple rank-5 right-hand side:
Gamma.U = [ exp(x).*sin(pi*x), exp(x).*sin(2*pi*x), exp(x).*sin(3*pi*x), exp(x).*sin(4*pi*x), exp(x).*sin(5*pi*x) ];
Gamma.S = diag([1,2,4,8,16]);
Gamma.V = [ exp(-2*y).*sin(pi*y), exp(-2*y).*sin(2*pi*y), exp(-2*y).*sin(3*pi*y), exp(-2*y).*sin(4*pi*y), exp(-2*y).*sin(5*pi*y) ];

end
