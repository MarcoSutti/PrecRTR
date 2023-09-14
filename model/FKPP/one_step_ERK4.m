function [ W1 ] = one_step_ERK4( W0, h, fun_RHS )

% function [ W1 ] = one_step_ERK4( W0, h, fun_RHS )
% Purpose: Perform one step of the explicit fourth-order Runge-Kutta (ERK4)
%          scheme.

% Created:     2023.05.09
% Last change: 2023.05.09

%   May 9, 2023:
%       Created.

K1 = fun_RHS( W0 );
K2 = fun_RHS( W0 + (h/2)*K1 );
K3 = fun_RHS( W0 + (h/2)*K2 );
K4 = fun_RHS( W0 + h*K3 );
W1 = W0 + (h/6)*( K1 + 2*K2 + 2*K3 + K4 );

end
