function [ F ] = cost_FKPP_f_from_struct( W, W0, W1, pars )

% function [ F ] = cost_FKPP_f_from_struct( W, W0, W1, pars )
% Purpose: Computes the cost function value for the Fisher-KPP equation
%          using the factorised form of W.
% Created:     2023.04.03
% Last change: 2023.04.28

%   Apr 28, 2023:
%       Use the factorized formats.
%       Checked coherence with the cost function in dense format.
%   Apr 27, 2023:
%       New cost function.
%   Apr 3, 2023:
%       Created.

%--------------------------------------------------------------------------
% Factorized Hadamard power of W1:
W1_sq = getFactorizedHadamard( W1, W1 );
%--------------------------------------------------------------------------

% 2023.04.28. Cost function with only the terms that depend on W:
G_W = W.U * W.S;
G_W0 = W0.U * W0.S;
G_W1 = W1.U * W1.S;
G_W1_sq = W1_sq.U * W1_sq.S;

sums_terms(1) = .5 * trace( (G_W'*pars.MmtMm)*G_W );
sums_terms(2) = -trace( (W.V'*W0.V) * (G_W0' * pars.MptMm * G_W) );
sums_terms(3) = -2*pars.dt * trace( (G_W1'*pars.Mminus * G_W) * (W.V' * pars.Romega * W1.V) );
sums_terms(4) = 2*pars.dt * trace( (G_W1_sq'*pars.Mminus * G_W) * (W.V' * pars.Romega *  W1_sq.V) );

F = sum(sums_terms);

end