function [ Uk, Sk, Vk ] = getTruncatedSVD( X, k )

% Returns the SVD of X truncated to rank k.
% Created:     20.02.2019
% Last change: 30.03.2020
%   March 26, 2020:
%       Replaced the version with the zero filling with a true truncated
%       version.

[ U, S, V ] = svd(X,0);

Uk = U(:,1:k);
Sk = S(1:k,1:k);
Vk = V(:,1:k);

end