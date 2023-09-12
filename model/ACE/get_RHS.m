function G = get_RHS( W, A, epsilon )

% function G = get_RHS( W, A, epsilon )

G = epsilon * (A * W + W * A) + W - W.^3;

end