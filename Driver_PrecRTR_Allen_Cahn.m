%==========================================================================
% Driver for using Preconditioned Riemannian Trust Region on the manifold
% of fixed-rank matrices.
% Script for ACE.
% Created:     2022.12.15
% Last change: 2023.02.03

%   Feb 3, 2023:
%       Cleanup.
%   Jan 30, 2023:
%       Introduced rank adaptivity.
%   Jan 23, 2023:
%       Corrected the formulation for L with periodic BCs.
%   Jan 10, 2023:
%       Added warm start with 'ACE_ref_sol_2023_W_t_1.mat' and the variables
%       t0, T, and Nt.
%   Dec 30, 2022:
%       Fixed mistake with periodic BCs.
%   Dec 15, 2022:
%       Created.
%==========================================================================
close all; clear; clc;

% Run startup_PrecRTR before.

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
pars.problem_type = 'ACE';

% Preconditiong: 0 for No, 1 for Yes.
pars.precon = 1;

options_lrie.verbosity = 1;
options_lrie.plot = false;
%--------------------------------------------------------------------------
% Diffusion coefficient of the Allen-Cahn equation:
pars.epsilon = 0.1;   % This is the same value used by Rodgers and Venturi.
%--------------------------------------------------------------------------
% Time step for numerical integration:
pars.dt = 0.20;
%--------------------------------------------------------------------------
% Size of the discretization
pars.lev_finest = 8;
pars.lev_coarsest = pars.lev_finest;
%--------------------------------------------------------------------------
% Manifold parameters
% Retraction type: 'metric', 'metric_explicit_inverse', 'ortho_ao', 'ortho_zhang'.
pars.retr_type = 'metric';

% Tolerance for rank truncation:
tol_rank = 1e-11;
%--------------------------------------------------------------------------
% Parameters for trustregions (Manopt)
options_tr.maxiter = 100;
% options_tr.minstepsize = 1e-14;
options_tr.tolgradnorm = 1e-9;
% options_tr.Delta_bar = 1;     % MS, 15.01.2021: It is better to comment
%                               out this line and let manopt take care of
%                               the Delta_bar automatically.
options_tr.verbosity = 0;
%--------------------------------------------------------------------------
% Spatial discretization:
Nx = 2^pars.lev_finest;
Ny = Nx;
Lx = 2*pi;
Ly = 2*pi;
hx = Lx/Nx;
hy = Ly/Ny;
x = linspace(-0.5*Lx+hx, 0.5*Lx, Nx+1); % !!!!!!!! Grazie JH !!!!!!!!
y = linspace(-0.5*Ly+hy, 0.5*Ly, Ny+1);
% x = linspace( 0, Lx, Nx);
% y = linspace( 0, Ly, Ny);
x = x(1:end-1); % !!!!!!!! Grazie JH !!!!!!!!
y = y(1:end-1);
[ xx , yy ] = ndgrid(x,y);
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf(['| Problem type: ', pars.problem_type, '                                            |\n'])
fprintf(['| Finest level: ', num2str(pars.lev_finest), ...
    '  (points: ', num2str(Nx), 'x', num2str(Nx), ')                           |\n'])
% fprintf(['| Low rank: ', num2str(pars.K), '                                                  |\n'])
if pars.precon==0
    fprintf('| WITHOUT PRECONDITIONER                                       |\n')
elseif pars.precon==1
    fprintf('| WITH PRECONDITIONER                                          |\n')
end
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                     Precompute Ah and L                      |\n');
fprintf('+--------------------------------------------------------------+\n');

global Precomputed;
Precomputed = struct( 'L', {}, 'Ah', {} );

%--------------------------------------------------------------------------
% Create matrix L in sparse format
% The matrix L is used to calculate forward finite differences of the
% first derivatives
B = [ -ones(Nx,1), [0; ones(Nx-1,1)] ];
d = [0,1];
Precomputed(1).L = (1/hx)*spdiags( B, d, Nx, Nx ); %(-eye(n_h)+diag(ones(n_h-1,1),1));

% Impose periodic BCs on L:
Precomputed(1).L(end,1) = Precomputed(1).L(1,2);

% % Get the second-order periodic spectral differentiation matrix:
% column = [ -pi^2/(3*hx^2)-1/6 -.5*(-1).^(1:Nx-1)./sin(hx*(1:Nx-1)/2).^2 ];
% D2 = (2*pi/Lx)^2 * toeplitz(column);   % 2022.12.30: avevo sbagliato il prefactor!!!

% The discretized minus Laplacian (without the prefactor 1/h^2) in sparse format:
Precomputed(1).Ah = get_Ah( Nx );
%--------------------------------------------------------------------------
% % Initial condition from the paper of Rodgers and Venturi, 2022:
% Wh0_mat = get_ACE_IC_RV2022( xx, yy );
% %--------------------------------------------------------------------------
% 
% % Truncate Wh0 to low rank:
% [ U, S, V ] = svd(Wh0_mat);
% 
% boolean_vect = diag(S) > tol_rank;
% 
% % Define the new rank:
% pars.K = nnz(boolean_vect);
% 
% % Truncate
% W.U = U(:,1:pars.K);
% W.S = S(1:pars.K,1:pars.K);
% W.V = V(:,1:pars.K);
% %--------------------------------------------------------------------------
% % Parameters for time evolution:
% t0 = 0;
% T = 15;
% Nt = ( T - t0 ) / pars.dt;
%--------------------------------------------------------------------------
% MS, 2023.01.10:
% Faccio partire la simulazione da t = 0.5 s, voglio vedere se cambia
% qualcosa nell'evoluzione temporale.
load("ACE_ref_256x256_T0.5_dt0.0001.mat")

[ U, S, V ] = svd( W_hist(:,:,end) );

boolean_vect = diag(S) > tol_rank;

% Define the new rank:
pars.K = nnz(boolean_vect);

% Truncate
W.U = U(:,1:pars.K);
W.S = S(1:pars.K,1:pars.K);
W.V = V(:,1:pars.K);
%--------------------------------------------------------------------------
% Parameters for time evolution:
t0 = t_hist(end);
T = 15.1;
Nt = ( T - t0 ) / pars.dt;
%--------------------------------------------------------------------------

% Pick the manifold of matrices of size n_h x n_h of fixed rank k.
problem.M = fixedrankembeddedfactory( Nx, Ny, pars.K, pars.retr_type );
problem.cost = @(X) cost_ACE_f_from_struct( X, W, pars );
problem.egrad = @(X) egrad_ACE_f_from_struct( X, W, pars );
problem.ehess = @(X,H) ehess_ACE_f_from_struct( X, H, pars );

%--------------------------------------------------------------------------
% Check gradient and Hessian:
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                         Check gradient                       |\n');
fprintf('+--------------------------------------------------------------+\n');
checkgradient(problem);
pause(1)
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                          Check Hessian                       |\n');
fprintf('+--------------------------------------------------------------+\n');
checkhessian(problem);
pause(1)
return
%--------------------------------------------------------------------------

fprintf('+--------------------------------------------------------------+\n');
fprintf('|       Start time integration of Allen-Cahn equation...       |\n');
fprintf('+--------------------------------------------------------------+\n');

% Initialization of t_hist_lra, rank_W_hist_lra, and W_hist_lra:
t_hist_lra = zeros( 1, Nt+1 );
rank_W_hist_lra = zeros( 1, Nt+1 );
W_hist_lra = struct( 'U', {}, 'S', {}, 'V', {} );

t_hist_lra(1) = t0;
W_hist_lra(1).U = W.U;
W_hist_lra(1).S = W.S;
W_hist_lra(1).V = W.V;
rank_W_hist_lra(1) = pars.K;

if options_lrie.verbosity >= 1
    fprintf('Rank of the initial condition: %d \n', pars.K );
end

for iter=1:Nt

    if options_lrie.verbosity >= 1
        fprintf('Iter: %.1d \n', iter);
    end

    % Define the new cost function and gradient:
    problem.M = fixedrankembeddedfactory( Nx, Ny, pars.K, pars.retr_type );
    problem.cost = @(X) cost_ACE_f_from_struct( X, W, pars );
    problem.egrad = @(X) egrad_ACE_f_from_struct( X, W, pars );
    problem.ehess = @(X,H) ehess_ACE_f_from_struct( X, H, pars );

    % MS, added 10.01.2021. Define a preconditioner for the Hessian.
    if pars.precon==1
        problem.precon = @(X,H) getXi( X, H, pars );
    end

    if options_tr.verbosity > 0
        fprintf('+--------------------------------------------------------------+\n');
        fprintf('|               Start Riemannian Trust Regions...              |\n');
        fprintf('+--------------------------------------------------------------+\n');
    end

    % Solve the problem with Riemannian Trust Regions
    [ W, ~, info_tr, ~ ] = trustregions( problem, W, options_tr );
 
    current_t = iter*pars.dt + t0;

    % Store history of t and W for later postprocessing:
    t_hist_lra(iter+1) = current_t;
    W_hist_lra(iter+1).U = W.U;
    W_hist_lra(iter+1).S = W.S;
    W_hist_lra(iter+1).V = W.V;

    %----------------------------------------------------------------------
    % Rank adaption:
    boolean_vect = diag(W.S) > tol_rank;

    % Define the new rank:
    pars.K = nnz(boolean_vect);

    % Truncate
    W.U = W.U(:,1:pars.K);
    W.S = W.S(1:pars.K,1:pars.K);
    W.V = W.V(:,1:pars.K);
    %----------------------------------------------------------------------

    rank_W_hist_lra(iter+1) = pars.K;

    % time_vect(irun) = toc;
    tr_gradnorm = [info_tr.gradnorm];
    num_inner_array = [info_tr.numinner];
    ntot_iter = length(tr_gradnorm) - 1;  % do not count the "0th" iteration
    % fprintf('Total time: %.2f s.\n', time_vect(irun) );

    if options_lrie.verbosity >= 2
        fprintf('Last grad. norm: %.4e.\n', tr_gradnorm(end) );
        fprintf('Residual: %.4e (norm of the Euclidean gradient).\n', norm_struct( problem.egrad( W ) ) );
        % CFR norm_Gh = norm_struct( getTruncationStruct( problem.egrad( Wh ), pars.K ) )
        % end

        fprintf('****************************************************************\n');
        fprintf('Number of outer iterations: %.2d\n', ntot_iter );
        fprintf('Sum number of inner iters:  %.2d\n', sum(num_inner_array) );
        fprintf('Max number of inner iters:  %.2d\n', max(num_inner_array) );
        fprintf('****************************************************************\n');
    end

    if options_lrie.plot
        % Form the full-rank matrix only for plotting purposes:
        W_mat = W.U * W.S * W.V';
        [ ~, handle_cf ] = contourf(x, y, real(W_mat), 100);
        set( handle_cf, 'edgecolor', 'none' );
        axis equal
        colormap("parula");
        colorbar;
        title(['Low-rank+IE, iter ',num2str(iter),', time = ', num2str(current_t)])
        drawnow;
    end

end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%--------------------------------------------------------------------------
% SAVE DATA TO MAT-FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save data                          |\n');
fprintf('+--------------------------------------------------------------+\n');

fileName_mfile = [ 'mdatafiles/', pars.problem_type, '_Prec', num2str(pars.precon), ...
    '_W_hist_lra_', num2str(Nx), 'x', num2str(Nx), '_K', num2str(pars.K) '_dt', ...
    num2str(pars.dt), '.mat'];

% Save all the three factors U, S, V separately:
save( fileName_mfile, 't_hist_lra', 'W_hist_lra', 'rank_W_hist_lra' );
fprintf('Saved data to file %s.\n', fileName_mfile);
%--------------------------------------------------------------------------