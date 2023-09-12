function [ M, rank_K ] = trunc_eps( M, eps_trunc )
% Use truncated SVD to truncate M to a low-rank matrix according to the
% tolerance eps_trunc.

[ U, S, V ] = svd(M);

boolean_vect = diag(S) > eps_trunc;

% Define the new rank:
rank_K = nnz(boolean_vect);

% Truncate M to rank_K:
M = U(:,1:rank_K) * S(1:rank_K,1:rank_K) * V(:,1:rank_K)';

end