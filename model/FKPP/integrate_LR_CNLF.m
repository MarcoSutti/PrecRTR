function [ Wn_RTR_hist, rank_W_RTR_hist ] = integrate_LR_CNLF( W0, Wn_minus_1, pars, options_RTR )

% function [ Wn_RTR_hist, rank_W_RTR_hist ] = integrate_LR_CNLF( W0, Wn_minus_1, pars, options_RTR )
% Purpose: Perform LR-CNLF.
%          Uses the preconditioned Riemannian trust-region method.

% Created:     2023.05.09
% Last change: 2023.05.09

%   May 9, 2023:
%       Created.

fprintf('+--------------------------------------------------------------+\n');
fprintf('|                  Riemannian Trust Regions...                 |\n');
fprintf('+--------------------------------------------------------------+\n');

% Truncate the initial condition:
[ U, S, V ] = svd( W0 );

boolean_vect = diag(S) > pars.tol_rank;

% Define the new rank:
pars.K = nnz(boolean_vect);

Wn_RTR_hist(1).U = U(:,1:pars.K);
Wn_RTR_hist(1).S = S(1:pars.K,1:pars.K);
Wn_RTR_hist(1).V = V(:,1:pars.K);

rank_W_RTR_hist(1) = pars.K;

% Truncate the first approximate solution (the one computed with ERK4):
[ U, S, V ] = svd( Wn_minus_1 );
Wn_RTR_hist(2).U = U(:,1:pars.K);
Wn_RTR_hist(2).S = S(1:pars.K,1:pars.K);
Wn_RTR_hist(2).V = V(:,1:pars.K);

rank_W_RTR_hist(2) = pars.K;

%--------------------------------------------------------------------------
% 2.3 Riemannian Trust-Region Method
%--------------------------------------------------------------------------
for time_iter=3:pars.Nt
    
    % Dummy variables:
    Wn_minus_1_RTR.U = Wn_RTR_hist(time_iter-1).U;
    Wn_minus_1_RTR.S = Wn_RTR_hist(time_iter-1).S;
    Wn_minus_1_RTR.V = Wn_RTR_hist(time_iter-1).V;

    Wn_minus_2_RTR.U = Wn_RTR_hist(time_iter-2).U;
    Wn_minus_2_RTR.S = Wn_RTR_hist(time_iter-2).S;
    Wn_minus_2_RTR.V = Wn_RTR_hist(time_iter-2).V;


    % Define the cost function and gradient:
    problem.M = fixedrankembeddedfactory( pars.Nx, pars.Nr, pars.K, pars.retr_type );
    problem.cost = @(X) cost_FKPP_f_from_struct( X, Wn_minus_2_RTR, Wn_minus_1_RTR, pars );
    problem.egrad = @(X) egrad_FKPP_f_from_struct( X, Wn_minus_2_RTR, Wn_minus_1_RTR, pars );
    problem.ehess = @(X,H) ehess_FKPP_f_from_struct( X, H, pars );

    % MS, added 2023.05.04. Define a preconditioner for the Hessian.
    if pars.precon==1
        problem.precon = @(X,H) getXi_FKPP_precon( X, H, pars );
    end

    %         %----------------------------------------------------------------------
    %         % Check gradient and Hessian:
    %         figure(1)
    %         fprintf('+--------------------------------------------------------------+\n');
    %         fprintf('|                         Check gradient                       |\n');
    %         fprintf('+--------------------------------------------------------------+\n');
    %         checkgradient(problem);
    %         pause(.5)
%     
    %     figure(2)
    %     fprintf('+--------------------------------------------------------------+\n');
    %     fprintf('|                          Check Hessian                       |\n');
    %     fprintf('+--------------------------------------------------------------+\n');
    %     checkhessian(problem);
    %     drawnow;
    %     %     pause(.5)
    %     pause
    %     %----------------------------------------------------------------------

    if options_RTR.verbosity > 0
        fprintf('+--------------------------------------------------------------+\n');
        fprintf('|               Start Riemannian Trust Regions...              |\n');
        fprintf('+--------------------------------------------------------------+\n');
        if pars.precon==0
            fprintf('| WITHOUT PRECONDITIONER                                       |\n')
        elseif pars.precon==1
            fprintf('| WITH PRECONDITIONER                                          |\n')
        end
        fprintf('+--------------------------------------------------------------+\n');
    end

    % Solve the problem with Euclidean Trust Regions
    [ Wn_RTR, ~, ~, ~ ] = trustregions( problem, Wn_minus_1_RTR, options_RTR );

    % Store the new approximate solution:
    current_t = (time_iter-1) * pars.dt;
    fprintf( "RTR: Time: %.4f \n", current_t );
%     t_hist(time_iter) = current_t;

%     Wn_RTR_hist(time_iter).U = Wn_RTR.U;
%     Wn_RTR_hist(time_iter).S = Wn_RTR.S;
%     Wn_RTR_hist(time_iter).V = Wn_RTR.V;

    %----------------------------------------------------------------------
    % Rank adaption:
    % Penso che dobbiamo tenere lo stesso rango ogni due, altrimenti non
    % posso fare certe operazioni.
    boolean_vect = diag(Wn_RTR.S) > pars.tol_rank;

    % Define the new rank:
    pars.K = nnz(boolean_vect);

    % Truncate
    Wn_RTR_hist(time_iter).U = Wn_RTR.U(:,1:pars.K);
    Wn_RTR_hist(time_iter).S = Wn_RTR.S(1:pars.K,1:pars.K);
    Wn_RTR_hist(time_iter).V = Wn_RTR.V(:,1:pars.K);
    %----------------------------------------------------------------------

    rank_W_RTR_hist(time_iter) = pars.K;

end