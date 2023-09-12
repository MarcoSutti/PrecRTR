function [ xi ] = getXi_FKPP( X, eta, pars )

% function [ xi ] = getXi_FKPP( X, eta, pars )
% Purpose: Recover the tangent vector xi from the equation:
%                       P_X * ( A * xi ) * P_X = eta,
%          without inverting the big matrix, for the manifold of fixed-rank
%          matrices.
% Created:     2023.04.07
% Last change: 2023.04.13

%   Apr 7, 2023:
%       Created by copying 'getXi.m'.

global Precomputed;

In = speye(pars.Nx);
Ik = speye(pars.K);

A_tilde = pars.hx * (-pars.dt * Precomputed(1).Ah + In);

% Get the spectral decomposition of X.U'*Ah*X.U
[ Q, D ] = eig( X.U'*A_tilde*X.U );
Q = Q';

% Transformations by Q
M_eta_hat = Q * eta.M;
U_hat = X.U * Q';
Vp_eta_hat = eta.Vp * Q';

% Frequently used quantities:
A_tilde_U_hat = A_tilde * U_hat;
U_hat_t_A_tilde = U_hat' * A_tilde;

% Orthogonal projectors % n-by-r
P_U_perp_A_tilde_U_hat = A_tilde_U_hat - X.U * (X.U'*A_tilde_U_hat);

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
    T_tilde_i = [ D(i,i) * In, X.V;
        X.V', zeros(pars.K,pars.K)];

    c1 = [ Vp_eta_hat(:,i); zeros(pars.K,1) ];

    T_tilde_inv_c1(:,i) = T_tilde_i\c1;

    Vp_xi_hat(:,i) = T_tilde_inv_c1(1:pars.Nx,i);

end

% Store in xi.Up:
xi.Up = Up_xi;

% Undo transformation done by Q and store in xi.Vp:
xi.Vp = Vp_xi_hat*Q;

end