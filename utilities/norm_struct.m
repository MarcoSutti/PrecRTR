function [ norm_X ] = norm_struct( X )

% Purpose: Computes Frobenius norm of X from its fatorized format, using
%          the trace trick.
% Created:     25.06.2019
% Last change: 20.04.2020

G = X.U * X.S;
H = X.V;

[ ~, Rg ] = qr(G,0);
[ ~, Rh ] = qr(H,0);

norm_X = sqrt( trace( (Rg'*Rg) * (Rh'*Rh) ) );

% norm_X = sqrt( trace( (G'*G) * (H'*H) ) );

end