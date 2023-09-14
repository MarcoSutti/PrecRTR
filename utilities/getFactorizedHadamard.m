function [ X_Hadamard_Y ] = getFactorizedHadamard( X, Y )

% Purpose: Given two X and Y matrices in the factorized form USV', computes
%          the Hadamard product X.*Y in the same factorized format, i.e.,
%          without forming the full matrices.
% Reference: Kressner and Tobler, 2014, "Algorithm 941: Htucker—a matlab
%            toolbox for tensors in hierarchical Tucker format", §7.
% Created:     19.06.2019
% Last change: 19.03.2020

%   March 19, 2020:
%       It seems that the Khatri-Rao product implemented by Laurent Sorber
%       (Laurent.Sorber@cs.kuleuven.be), called kr.m, is more efficient
%       than the version implemented in the MATLAB Tensor Toolbox (called
%       khatrirao.m).

% Use the kr.m implemented by Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
X_Hadamard_Y.U = kr( X.U', Y.U' )';
X_Hadamard_Y.S = kron( X.S, Y.S );
X_Hadamard_Y.V = kr( X.V', Y.V')';

% Previous version used the khatrirao.m function from the MATLAB Tensor
% Toolbox
% X_Hadamard_Y.U = khatrirao( X.U', Y.U' )';
% X_Hadamard_Y.S = kron( X.S, Y.S );
% X_Hadamard_Y.V = khatrirao( X.V', Y.V')';

end