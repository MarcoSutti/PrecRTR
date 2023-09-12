function [ M_r ] = trunc_r( M, r )
% Use truncated SVD to truncate M to a low-rank matrix according to the
% fixed rank r

[ U, S, V ] = svd(M);

% Truncate M to rank r and form again the matrix
M_r = U(:,1:r) * S(1:r,1:r) * V(:,1:r)';

end