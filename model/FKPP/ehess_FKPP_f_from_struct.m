function [ EH ] = ehess_FKPP_f_from_struct( W, H, pars )

% function [ EH ] = ehess_FKPP_f_from_struct( W, H, pars )
% Purpose: Computes the directional derivative of the gradient along the
%          direction of H for the Fisher-KPP equation.
% Created:     2023.04.03
% Last change: 2023.04.28

%   Apr 28, 2023:
%       Use the factorized formats.
%       Checked coherence with the Euclidean Hessian in dense format.
%   Apr 27, 2023:
%       New cost function.
%   Apr 3, 2023:
%       Created.

%--------------------------------------------------------------------------
% Convert H to the ambient space format:
H_amb.U = [W.U*H.M + H.Up, W.U];
[ ~, k_h ] = size(H_amb.U);
H_amb.S = speye( k_h );
H_amb.V = [W.V, H.Vp];

% Computes the Euclidean Hessian in factorized format:
EH.U = pars.MmtMm * H_amb.U;
EH.S = H_amb.S;
EH.V = H_amb.V;

end