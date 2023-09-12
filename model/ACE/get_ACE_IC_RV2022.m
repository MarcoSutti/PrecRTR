function [ Wh0_mat ] = get_ACE_IC_RV2022( xx, yy )

% function [ Wh0_mat ] = get_ACE_IC_RV2022( xx, yy )
% Purpose: Calculate the second-order periodic spectral differentiation
% matrix of [Tre00, eq. (3.12)] rescaled to [−L, L].

% Reference: [RV22] Rodgers, A. and Venturi, D., 2022.

% Created:     2023.02.03
% Last change: 2023.03.03

%   Feb 3, 2023:
%       Created.

% Initial condition for the  from the paper of Rodgers and Venturi, 2022:
Wh0_mat = get_U(xx, yy) - get_U(xx, 2*yy) + get_U(3*xx+pi, 3*yy+pi) ...
    - 2*get_U(4*xx, 4*yy) + 2*get_U(5*xx, 5*yy);

end


function U = get_U( X, Y )

% Compute the inital condition f_0(x,y) according to equation (77):
% Computes u(x,y) according to eq. (78).

num = ( exp(-tan(X).^2) + exp(-tan(Y).^2) ) .* sin(X) .* sin(Y);
denom = 1 + exp(abs(csc(-X./2))) + exp(abs(csc(-Y./2)));

U = num ./ denom;

end