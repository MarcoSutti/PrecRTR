function [ F ] = cost_ACE_f_from_struct( W, W_tilde, pars )

% function [ F ] = cost_ACE_f_from_struct( W, W_tilde, pars )
% Purpose: Computes the cost function value for the Allen-Cahn equation
%          using the factorised form of W.
% Created:     2022.12.19
% Last change: 2023.01.09

%   Dec 19, 2022:
%       Created.

global Precomputed;

W_mat = W.U * W.S * W.V';
% W_tilde_mat = W_tilde.U * W_tilde.S * W_tilde.V';

G_1 = W.U * W.S;
H_1 = W.V;

G_2 = W_tilde.U * W_tilde.S;
H_2 = W_tilde.V;
%--------------------------------------------------------------------------
Nx = size(W.U,1);
Lx = 2*pi;
% Ly = 2*pi;
hx = Lx/Nx;
h_sq = hx^2;
%--------------------------------------------------------------------------
% Compute the factorized truncated elementwise product Wh.*Wh
% 04.12.2019: Keep the Wh_sq of the size k^2-by-k^2, i.e., do not truncate it!
% Wh_sq = getFactorizedHadamard( Wh, Wh );
% 
% G_2 = Wh_sq.U * Wh_sq.S;
% H_2 = Wh_sq.V;

% 1) Get the first-order periodic spectral differentiation matrix:
% column = [0 .5*(-1).^(1:Nx-1).*cot((1:Nx-1)*hx/2)];
% D1 = toeplitz(column,column([1 Nx:-1:2]));

% column = [0 .5*(-1).^(1:N-1).*cot((1:N-1)*h/2)]';
% D = toeplitz(column,column([1 N:-1:2]));

% Factorized first derivatives
LG = Precomputed(1).L*G_1;
% LG = D1 * G_1;

LH = Precomputed(1).L*H_1;
% LH = D1 * H_1;

GtG = G_1'*G_1;
HtH = H_1'*H_1;

pref = .5 * pars.epsilon * pars.dt;

first_term = pref * trace( (LG'*LG) * HtH + GtG * (LH'*LH) );

second_term = .5 * (1-pars.dt) * sum(diag(W.S).^2);

% Terzo termine: Per il momento lo lascio così:
third_term = .25 * pars.dt * sum(sum(W_mat.^4));

fourth_term = - trace( (G_2'*G_1) * (H_1'*H_2) );

F = h_sq * ( first_term + second_term + third_term + fourth_term );


% forward difference
% F = 0.5 * area *  ...

% % Check:
% first_term_mat = pref * ( norm(D1 * W_mat, 'fro')^2 + norm(W_mat * D1', 'fro')^2 );
% 
% second_term_mat = (1-pars.dt)/2 * sum(sum(W_mat.^2));
% 
% third_term_mat = pars.dt/4 * sum(sum(W_mat.^4));
% 
% fourth_term_mat = -sum(sum(W_tilde_mat .* W_mat ));
% 
% % 
% norm( first_term  - first_term_mat, 'fro')

% % norm( second_term - second_term_mat, 'fro')
% % norm( third_term  - third_term_mat, 'fro')
% norm( fourth_term - fourth_term_mat, 'fro' )

end