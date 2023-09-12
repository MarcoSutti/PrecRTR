function [ xi ] = getXi_FKPP_precon( W, eta, pars )

% function [ xi ] = getXi_FKPP_precon( W, eta, pars )
% Purpose: Recover the tangent vector xi from the equation:
%                       P_X * ( A * xi ) * P_X = eta,
%          without inverting the big matrix, for the manifold of fixed-rank
%          matrices.
% Created:     2023.04.18
% Last change: 2023.05.08


%   May 6, 2023:
%       A_tilde should be MmtMm. See my latexed notes on "Fisher–KPP PDE".
%   Apr 18, 2023:
%       Created.

In = speye(pars.Nx);
Ik = speye(pars.K);

A_tilde = pars.MmtMm;

% Get the spectral decomposition of X.U'*Ah*X.U
[ Q, D ] = eig( W.U'*A_tilde*W.U );
Q = Q';

% Transformations by Q
M_eta_hat = Q * eta.M;
U_hat = W.U * Q';
Vp_eta_hat = eta.Vp * Q';

% Frequently used quantities:
A_tilde_U_hat = A_tilde * U_hat;
U_hat_t_A_tilde = U_hat' * A_tilde;

% Orthogonal projectors % n-by-r
P_U_perp_A_tilde_U_hat = A_tilde_U_hat - W.U * (W.U'*A_tilde_U_hat);

% Initialize:
vi = zeros( pars.K );
K = cell(1,pars.K);

% They do not depend on i in this case:
XiU_2nd_term = A_tilde\U_hat;
SiU = U_hat' * XiU_2nd_term;

for i=1:pars.K
    %----------------------------------------------------------------------
    % Find the vi's and K{i}'s
    %----------------------------------------------------------------------

    BiU = [ eta.Up(:,i),  P_U_perp_A_tilde_U_hat ];
    
    XiU_1st_term = A_tilde\BiU;
    
    NiU = SiU\(U_hat' * XiU_1st_term);
    
    XiU = XiU_1st_term - XiU_2nd_term * NiU;

    vi(:,i) = U_hat_t_A_tilde * XiU(:,1);
    
    K{i} = - U_hat_t_A_tilde * XiU(:,2:end);
end

% Form the block-diagonal matrices calligraphic K and K_tilde
% Form the block-diagonal matrices calligraphic K and K_tilde
calK = sparse( blkdiag( K{:} ) );

% Form the big matrix of the linear system for M_xi_hat
bigMat = calK + sparse( kron(Ik, D) );

% Form the RHS matrix R
R = M_eta_hat - vi;

vecM_xi_hat = bigMat\R(:);

M_xi_hat = reshape( vecM_xi_hat, [pars.K, pars.K]);

% Undo the transformations by Q and Q_tilde to recover xi.M
xi.M = Q'*M_xi_hat;

% Initialize Up_xi
Up_xi = zeros( pars.Nx, pars.K );
Vp_xi_hat = zeros( pars.Nx, pars.K );
T_tilde_inv_c1 = zeros( pars.Nx + pars.K, pars.K );

%--------------------------------------------------------------------------
% Recover Up_xi_hat and Vp_xi_hat
%--------------------------------------------------------------------------
for i=1:pars.K
    %----------------------------------------------------------------------
    % Recover Up_xi:
    %----------------------------------------------------------------------    
    biU = [ eta.Up(:,i),  P_U_perp_A_tilde_U_hat * M_xi_hat(:,i) ];
    
    xiU_1st_term = A_tilde\biU;
    
    % Use Cholesky factorization to solve the system with Si
    niU = SiU\(U_hat' * xiU_1st_term);
    
    XiU = xiU_1st_term - XiU_2nd_term * niU;
    
    Up_xi(:,i) = XiU(:,1) - XiU(:,2);

    %----------------------------------------------------------------------
    % Recover Vp_xi_hat:
    %----------------------------------------------------------------------
    T_tilde_i = [ D(i,i) * In, W.V;
        W.V', zeros(pars.K,pars.K)];

    c1 = [ Vp_eta_hat(:,i); zeros(pars.K,1) ];

    T_tilde_inv_c1(:,i) = T_tilde_i\c1;

    Vp_xi_hat(:,i) = T_tilde_inv_c1(1:pars.Nx,i);

end

% Store in xi.Up:
xi.Up = Up_xi;

% Undo transformation done by Q and store in xi.Vp:
xi.Vp = Vp_xi_hat*Q;

% % Check correctness of the preconditioner:
% eta_ambient = W.U * eta.M * W.V' + eta.Up * W.V' + W.U * eta.Vp';
% xi_ambient = W.U * xi.M * W.V' + xi.Up * W.V' + W.U * xi.Vp';
% 
% EH = pars.MmtMm * xi_ambient;
% 
% PU = W.U*W.U';
% PV = W.V*W.V';
% PUperp = eye(pars.Nx) - PU;
% PVperp = eye(pars.Nx) - PV;
% 
% % Projection of EH onto the tangent space to Mr at W:
% PW_A_tilde_xi = PU * EH * PV + PUperp * EH * PV + PU * EH * PVperp;
% % H_W_times_xi = PU * (pars.MmtMm * xi_ambient) * PV + PUperp * (pars.MmtMm * xi_ambient) * PV + PU * (pars.MmtMm * xi_ambient) * PVperp;
% 
% norm(PW_A_tilde_xi - eta_ambient, "fro")
% 
% % System (6.3)
% norm(W.U'*pars.MmtMm*W.U*xi.M + W.U'*pars.MmtMm*xi.Up - eta.M, "fro")
% norm(PUperp * pars.MmtMm * W.U * xi.M + PUperp * pars.MmtMm * xi.Up - eta.Up, "fro")
% norm(W.U'*pars.MmtMm*W.U*xi.Vp' - eta.Vp', "fro")
% 
% pause


end