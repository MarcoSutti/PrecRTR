function [ xi ] = getXi( X, eta, pars )

% function [ xi ] = getXi( X, eta, pars )
% Purpose: Recover the tangent vector xi from the equation:
%                       P_X * ( A * xi + xi * A ) * P_X = eta,
%          without inverting the big matrix, for the manifold of fixed-rank
%          matrices.
% Created:     2021.01.09
% Last change: 2023.01.11

%   Jan 11, 2023:
%       Replaced the switch with sparse(blkdiag(K{:})).
%   Sept 21, 2022:
%       Replaced the if statement with the switch, case, otherwise.
%   Jan 12, 2021:
%       Pass pars as input.
%       Added idx.
%   Jan 12, 2021:
%       Avoid performing the same calculations twice.
%       Use Cholesky factorization to solve the system with Si.
%   Jan 11, 2021:
%       Solve the saddle point system with the Schur complement trick.
%   Jan 9, 2021:
%       Created.

global Precomputed;

%--------------------------------------------------------------------------
n = size(X.U,1);
k = pars.K;
%--------------------------------------------------------------------------
% lev = log2(n-1);
% idx = lev - pars.lev_coarsest + 1;
%--------------------------------------------------------------------------

A = Precomputed(1).Ah;

%--------------------------------------------------------------------------

In = speye(n);
Ik = speye(k);

% Get the spectral decomposition of X.U'*Ah*X.U
[ Q, D ] = eig( X.U'*A*X.U );
Q = Q';

% disp('lo sto usando')

% Get the spectral decomposition of X.V'*Ah*X.V
[ Q_tilde, D_tilde ] = eig( X.V'*A*X.V );
Q_tilde = Q_tilde';

% Transformations by Q
M_eta_hat = Q * eta.M * Q_tilde';
Up_eta_hat = eta.Up * Q_tilde';
Vp_eta_hat = eta.Vp * Q';
U_hat = X.U * Q';
V_hat = X.V * Q_tilde';

A_U_hat = A * U_hat;
A_V_hat = A * V_hat;

U_hat_t_A = U_hat' * A;
V_hat_t_A = V_hat' * A;

% Orthogonal projectors
P_U_perp_A_U_hat = A_U_hat - X.U * (X.U'*A_U_hat);
P_V_perp_A_V_hat = A_V_hat - X.V * (X.V'*A_V_hat);


% Initialize:
vi = zeros( k );

vi_tilde = zeros( k );
K = cell(1,k);
K_tilde = cell(1,k);

A_plus_d_tilde = cell(1,k);
A_plus_d = cell(1,k);
XiU_2nd_term = cell(1,k);
XiV_2nd_term = cell(1,k);
R_SiU = cell(1,k);
R_SiV = cell(1,k);

for i=1:k
    %----------------------------------------------------------------------
    % Find the vi's and K{i}'s
    %----------------------------------------------------------------------
    A_plus_d_tilde{i} = A + D_tilde(i,i) * In;
    
    %     Ti = [ A_plus_di_tilde,  U_hat;  U_hat',   zeros(k) ];
    
    %     bi_1st_term = [ Up_eta_hat(:,i); zeros(k,1) ];
    
    %     sol_1st_term = Ti\bi_1st_term;
    
    %     vi(:,i) = U_hat_t_A * sol_1st_term(1:n)
    
    XiU_2nd_term{i} = A_plus_d_tilde{i}\U_hat;
    
    SiU = U_hat' * XiU_2nd_term{i};
    
    BiU = [ Up_eta_hat(:,i),  P_U_perp_A_U_hat ];
    
    XiU_1st_term = A_plus_d_tilde{i}\BiU;
    
    % Use Cholesky factorization to solve the system with Si
    R_SiU{i} = chol(real(SiU));
    
    NiU = R_SiU{i}\(R_SiU{i}'\(U_hat' * XiU_1st_term));
    
    XiU = XiU_1st_term - XiU_2nd_term{i} * NiU;

    vi(:,i) = U_hat_t_A * XiU(:,1);
    
    K{i} = D_tilde(i,i) * Ik - U_hat_t_A * XiU(:,2:end);
    
    %----------------------------------------------------------------------
    % Find the vi_tilde's and K_tilde{i}'s
    %----------------------------------------------------------------------
    A_plus_d{i} = A + D(i,i) * In;
    
    XiV_2nd_term{i} = A_plus_d{i}\V_hat;
    
    SiV = V_hat' * XiV_2nd_term{i};
    
    BiV = [ Vp_eta_hat(:,i),  P_V_perp_A_V_hat ];
    
    XiV_1st_term = A_plus_d{i}\BiV;
    
    % Use Cholesky factorization to solve the system with Si
    R_SiV{i} = chol(real(SiV));
    
    NiV = R_SiV{i}\(R_SiV{i}'\(V_hat' * XiV_1st_term));
    
    XiV = XiV_1st_term - XiV_2nd_term{i} * NiV;
    
    vi_tilde(:,i) = V_hat_t_A * XiV(:,1);
    
    K_tilde{i} = D(i,i) * Ik - V_hat_t_A * XiV(:,2:end);
    
end

% Form the block-diagonal matrices calligraphic K and K_tilde
calK = sparse( blkdiag( K{:} ) );
calK_tilde = sparse( blkdiag( K_tilde{:} ) );

Pkk = perfect_shuffle( k, k );

% Form the big matrix of the linear system for M_xi_hat
bigMat = calK + Pkk * calK_tilde * Pkk;

% Form the RHS matrix R
R = M_eta_hat - vi - vi_tilde';

vecM_xi_hat = bigMat\R(:);

M_xi_hat = reshape( vecM_xi_hat, [k, k]);

% Undo the transformations by Q and Q_tilde to recover xi.M
xi.M = Q'*M_xi_hat*Q_tilde;


% Initialize Up_xi_hat
Up_xi_hat = zeros( n, k );

% Initialize Vp_xi_hat
Vp_xi_hat = zeros( n, k );

%--------------------------------------------------------------------------
% Recover Up_xi_hat and Vp_xi_hat
%--------------------------------------------------------------------------
for i=1:k
    %----------------------------------------------------------------------
    % Recover Up_xi_hat:
    %----------------------------------------------------------------------    
    biU = [ Up_eta_hat(:,i),  P_U_perp_A_U_hat * M_xi_hat(:,i) ];
    
    xiU_1st_term = A_plus_d_tilde{i}\biU;
    
    % Use Cholesky factorization to solve the system with Si
    niU = R_SiU{i}\(R_SiU{i}'\(U_hat' * xiU_1st_term));
    
    XiU = xiU_1st_term - XiU_2nd_term{i} * niU;
    
    Up_xi_hat(:,i) = XiU(:,1) - XiU(:,2);
    
    %----------------------------------------------------------------------
    % Recover Vp_xi_hat
    %----------------------------------------------------------------------
    biV = [ Vp_eta_hat(:,i),  P_V_perp_A_V_hat * M_xi_hat(i,:)' ];
    
    xiV_1st_term = A_plus_d{i}\biV;
    
    % Use Cholesky factorization to solve the system with Si
    niV = R_SiV{i}\(R_SiV{i}'\(V_hat' * xiV_1st_term));
    
    XiV = xiV_1st_term - XiV_2nd_term{i} * niV;
    
    Vp_xi_hat(:,i) = XiV(:,1) - XiV(:,2);
    
end

% Undo the transformation by Q_tilde to recover xi.Up
xi.Up = Up_xi_hat * Q_tilde;

% Undo the transformation by Q to recover xi.Vp
xi.Vp = Vp_xi_hat * Q;

end