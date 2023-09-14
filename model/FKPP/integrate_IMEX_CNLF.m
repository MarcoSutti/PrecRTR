function [ t_EuTR_hist, W_EuTR_hist, t_hist, W_CNLF_hist, rank_W_CNLF_hist ] = integrate_IMEX_CNLF( W0, Wn_minus_1, pars, options_EuTR )

% function [ t_hist, W_CNLF_hist, rank_W_CNLF_hist ] = integrate_IMEX_CNLF( W0, Wn_minus_1, pars)
% Purpose: Perform numerical time integration with the IMEX-CNLF scheme.

% Created:     2023.05.09
% Last change: 2023.09.14

%   May 9, 2023:
%       Created.

fprintf('+--------------------------------------------------------------+\n');
fprintf('|           IMEX-CNLF method and Euclidean TR method           |\n');
fprintf('+--------------------------------------------------------------+\n');

% MS, 2023.05.09: We need to avoid to store the entire history.
% Save only a subset of the computed data. Use a stride.
stride = floor(pars.Nt/100);

Wn_minus_1_CNLF = Wn_minus_1;
Wn_minus_2_CNLF = W0;

idx_stride = 1;
t_hist(idx_stride) = 0;
W_CNLF_hist(:,:,idx_stride) = W0;
rank_W_CNLF_hist(idx_stride) = rank( W0 );

% For the Euclidean Trust Regions:
W_EuTR_hist = zeros(pars.Nx,pars.Nr,pars.Nt);
W_EuTR_hist(:,:,1) = W0;
W_EuTR_hist(:,:,2) = Wn_minus_1;
t_EuTR_hist(1) = 0;
t_EuTR_hist(2) = pars.dt;

for time_iter=3:pars.Nt
    %----------------------------------------------------------------------
    % 2.1. IMEX-CNLF method
    %      NB: For all the realizations at once.
    %----------------------------------------------------------------------
    % Dummy variables:
    %     Wn_minus_1_CNLF = W_CNLF_hist(:,:,time_iter-1);
    %     Wn_minus_2_CNLF = W_CNLF_hist(:,:,time_iter-2);

    % IMEX-CNLF method
    B = pars.Mplus*Wn_minus_2_CNLF + 2*pars.dt*Wn_minus_1_CNLF.*(1-Wn_minus_1_CNLF)*pars.Romega;
    Wn_CNLF = pars.Mminus\B;

    % Store the new approximate solution:
    %     current_t = (time_iter-1) * pars.dt;
    %         fprintf( "Time: %.4f \n", current_t );
    %     t_hist(time_iter) = current_t;
    %     W_CNLF_hist(:,:,time_iter) = Wn_CNLF;
    %     rank_W_CNLF_hist(time_iter) = rank( Wn_CNLF );

    %----------------------------------------------------------------------
    % 2.2 Euclidean Trust-Region Method
    %----------------------------------------------------------------------
    % Dummy variables:
    Wn_minus_1_EuTR = W_EuTR_hist(:,:,time_iter-1);
    Wn_minus_2_EuTR = W_EuTR_hist(:,:,time_iter-2);

    % Define the cost function and gradient:
    problem.M = euclideanfactory(pars.Nx, pars.Nr);
    problem.cost = @(X) cost_FKPP_CNLF( X, Wn_minus_2_EuTR, Wn_minus_1_EuTR, pars );
    problem.egrad = @(X) egrad_FKPP_CNLF( X, Wn_minus_2_EuTR, Wn_minus_1_EuTR, pars );
    problem.ehess = @(X,H) ehess_FKPP_CNLF( X, H, pars );

    %----------------------------------------------------------------------
    % Check gradient and Hessian:
    %     figure(1)
    %     fprintf('+--------------------------------------------------------------+\n');
    %     fprintf('|                         Check gradient                       |\n');
    %     fprintf('+--------------------------------------------------------------+\n');
    %     checkgradient(problem);
    % %     drawnow;
    %     pause(.5)
    %
    %     figure(2)
    %     fprintf('+--------------------------------------------------------------+\n');
    %     fprintf('|                          Check Hessian                       |\n');
    %     fprintf('+--------------------------------------------------------------+\n');
    %     checkhessian(problem);
    %     drawnow;
    %     pause
    %     %----------------------------------------------------------------------

    if options_EuTR.verbosity > 0
        fprintf('+--------------------------------------------------------------+\n');
        fprintf('|                Start Euclidean Trust Regions...              |\n');
        fprintf('+--------------------------------------------------------------+\n');
    end

    % Solve the problem with Euclidean Trust Regions
    [ Wn_EuTR, ~, ~, ~ ] = trustregions( problem, Wn_minus_1_EuTR, options_EuTR );

    % Store the new approximate solution:
    current_t = (time_iter-1) * pars.dt;

    t_EuTR_hist(time_iter) = current_t;
    W_EuTR_hist(:,:,time_iter) = Wn_EuTR;

    % MS, 2023.05.09: Save the numerical solution only every 'stride' steps.
    if mod(time_iter-1,stride)==0
        fprintf( "IMEX-CNLF: Time: %.4f \n", current_t );
        idx_stride = idx_stride + 1;
        t_hist(idx_stride) = current_t;
        W_CNLF_hist(:,:,idx_stride) = Wn_CNLF;
        rank_W_CNLF_hist(idx_stride) = rank( Wn_CNLF );
    end

    % Update!!!!!!!!!!!
    Wn_minus_2_CNLF = Wn_minus_1_CNLF;
    Wn_minus_1_CNLF = Wn_CNLF;

    %----------------------------------------------------------------------
    % Check that the solution given by the IMEX-CNLF scheme and the
    % solution given by the Euclidean TR method are the same:
    check = norm(Wn_EuTR - Wn_CNLF);
    fprintf('Check norm(wn_Eucl_tr - wn_CNLF): %.5e.\n', check );

    if check > 100*options_EuTR.tolgradnorm
        error('They are not solving the same problem!!!');
    end

end
end