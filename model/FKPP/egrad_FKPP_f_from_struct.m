function [ EG ] = egrad_FKPP_f_from_struct( W, W0, W1, pars )

% function [ EG ] = egrad_FKPP_f_from_struct( W, W0, W1, pars )
% Purpose: Computes the Euclidean gradient in factorized low-rank format
%          for the Fisher-KPP equation.
% Created:     2023.04.03
% Last change: 2023.04.28

%   Apr 28, 2023:
%       Use the factorized formats.
%       Checked coherence with the Euclidean gradient in dense format.
%   Apr 27, 2023:
%       New cost function.
%   Apr 3, 2023:
%       Created.

%--------------------------------------------------------------------------
% Factorized Hadamard power of W1:
W1_sq = getFactorizedHadamard( W1, W1 );
%--------------------------------------------------------------------------

% Computes the Euclidean gradient in factorized format:
EG.U = [ pars.MmtMm * W.U, pars.MptMm' * W0.U, pars.Mminus'*W1_sq.U, pars.Mminus'*W1.U ];
EG.S = blkdiag( W.S, -W0.S, 2*pars.dt*W1_sq.S, -2*pars.dt*W1.S );
EG.V = [ W.V, W0.V, pars.Romega*W1_sq.V, pars.Romega*W1.V ];

end