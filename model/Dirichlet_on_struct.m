function [ Wh ] = Dirichlet_on_struct( Wh )

% function [ Wh ] = Dirichlet_on_struct( Wh )
% Purpose: Apply homogeneous Dirichlet BCs directly on the factored form of
%          Wh, without forming the full matrix.
% Created:     08.04.2019
% Last change: 08.04.2019

Wh.U(1,:) = 0;
Wh.U(end,:) = 0;
Wh.V(1,:) = 0;
Wh.V(end,:) = 0;

end