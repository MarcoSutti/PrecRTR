%==========================================================================
% Driver for using Preconditioned Riemannian Trust Region on the manifold
% of fixed-rank matrices.
% Script for the 'LYAP' and the 'NPDE' problems.
% Created:     2020.03.02
% Last change: 2023.09.13

%   Sep 13, 2023:
%       Added option to switch on/off the checks on consistency of the
%       gradient and the Hessian.
%   Feb 20, 2023:
%       Added options_tr.maxinner for the inner tCG solver.
%   Sep 4, 2022:
%       Little change to plotting functions.
%   Aug 29, 2022:
%       Added saving of simulation data to mfiles.
%       Added startup.
%   Aug 25, 2022:
%       Added plotting utilities.
%   Jan 10, 2021:
%       Added Hessian preconditioner.
%   Dec 18, 2020:
%       Changed 'startup' to 'close all; clear; clc;' and 'RHS_LYAP' to
%       'RHS_LYAP_struct'
%   April 19, 2020:
%       Added loop for performing several runs.
%   April 6, 2020:
%       Changed the scaling factor in the initial guess, now it depends on
%       the discretization parameter h.
%   March 20, 2020:
%       Changed the way the initial guess is provided.
%==========================================================================
close all; clear; clc;

% Run startup_PrecRTR before.

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Fix stream of rng
rng(1)

% Number of runs to perform several times the same experiment:
number_of_runs = 1;
%--------------------------------------------------------------------------
% Problem type: 'NPDE', 'LYAP'.
pars.problem_type = 'LYAP';
%--------------------------------------------------------------------------
% Check the consistency of the gradient and the Hessian matrices with the 
% cost function.
% 0 for No, 1 for Yes.
check_grad_hess = 0;
%--------------------------------------------------------------------------
% Preconditiong: 0 for No, 1 for Yes.
pars.precon = 1;
%--------------------------------------------------------------------------
pars.lambda = 10;     % Only used if 'NPDE' is chosen as pars.problem_type.
%--------------------------------------------------------------------------
% Size of the discretization
pars.lev_finest = 15;
pars.lev_coarsest = pars.lev_finest;
%--------------------------------------------------------------------------
% Manifold parameters
% Retraction type: 'metric', 'metric_explicit_inverse', 'ortho_ao', 'ortho_zhang'.
pars.retr_type = 'metric';

% Set the low-rank
pars.K = 5;
%--------------------------------------------------------------------------
% Parameters for trustregions (Manopt)
options_tr.maxiter = 300;
% options_tr.maxinner = 1000;   % max number of iterations for the inner tCG solver
options_tr.minstepsize = 1e-12;
options_tr.tolgradnorm = 1e-12;
% options_tr.Delta_bar = 0.5;     % MS, 15.01.2021: It is better to comment
%                               out this line and let manopt take care of
%                               the Delta_bar automatically.
options_tr.verbosity = 0;
%--------------------------------------------------------------------------
% % For logging the command window output:
% fileName_log = [ 'logs/', pars.problem_type, '_Prec', num2str(pars.precon), ...
%     '_Grad_', num2str(pars.lev_finest), '_K', num2str(pars.K), '.log' ];
% 
% diary(fileName_log)
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf(['| Problem type: ', pars.problem_type, '                                           |\n'])
fprintf( '| Finest level: %d \n', pars.lev_finest )
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                Precompute Gamma_h, Ah and L                  |\n');
fprintf('+--------------------------------------------------------------+\n');

tic

global Gamma_h;
Gamma_h = struct( 'U', {}, 'S', {}, 'V', {} );

global Precomputed;
Precomputed = struct( 'L', {}, 'Ah', {} );

n_h = 2^pars.lev_finest + 1;

h = 1/(n_h-1);
%--------------------------------------------------------------------------
% Create matrix L in sparse format
% The matrix L is used to calculate forward finite differences of the
% first derivatives
B = [ -ones(n_h,1), [0; ones(n_h-1,1)] ];
d = [0,1];
Precomputed(1).L = (1/h)*spdiags( B, d, n_h, n_h ); %(-eye(n_h)+diag(ones(n_h-1,1),1));

% The discretized minus Laplacian (without the prefactor 1/h^2) in sparse format:
Precomputed(1).Ah = spdiags( ones(n_h,1) * [-1 2 -1], -1:1, n_h, n_h );
%--------------------------------------------------------------------------

% Impose homogeneous Dirichlet BCs on the discretized Laplacian:
Precomputed(1).Ah(1,:)=0;
Precomputed(1).Ah(end,:)=0;
Precomputed(1).Ah(:,1)=0;
Precomputed(1).Ah(:,end)=0;

%--------------------------------------------------------------------------
if strcmp(pars.problem_type, 'LYAP')
    Gamma_h_lev_i = RHS_LYAP_struct( h );
elseif strcmp(pars.problem_type, 'NPDE')
    Gamma_h_lev_i = RHS_NPDE_struct( h );
end
Gamma_h(1).U = Gamma_h_lev_i.U;
Gamma_h(1).S = Gamma_h_lev_i.S;
Gamma_h(1).V = Gamma_h_lev_i.V;

toc
%--------------------------------------------------------------------------
pars.lev = pars.lev_finest;
%--------------------------------------------------------------------------
% Pick the manifold of matrices of size n_h x n_h of fixed rank k.
problem.M = fixedrankembeddedfactory( n_h, n_h, pars.K, pars.retr_type );
if strcmp(pars.problem_type, 'LYAP')
    problem.cost = @(X) cost_LYAP_f_from_struct( X, pars );
    problem.egrad = @(X) egrad_LYAP_f_from_struct( X, pars );
    problem.ehess = @(X,H) ehess_LYAP_f_from_struct( X, H, pars );
elseif strcmp(pars.problem_type, 'NPDE')
    problem.cost = @(X) cost_NPDE_f_from_struct( X, pars );
    problem.egrad = @(X) egrad_NPDE_f_from_struct( X, pars );
    problem.ehess = @(X,H) ehess_NPDE_f_from_struct( X, H, pars );
end

% MS, added 10.01.2021. Define a preconditioner for the Hessian.
if pars.precon==1
    problem.precon = @(X,H) getXi( X, H, pars );
end

if check_grad_hess==1
    % Check gradient and Hessian:
    checkgradient(problem);
    pause(0.5)
    checkhessian(problem);
    pause(0.5)
end

fprintf('+--------------------------------------------------------------+\n');
fprintf('|                    Set initial point Wh0                     |\n');
fprintf('+--------------------------------------------------------------+\n');

tic

scaling_factor = h;

Wh0.U = scaling_factor * rand( n_h, pars.K );
Wh0.S = h^2 * eye( pars.K );
Wh0.V = scaling_factor * rand( n_h, pars.K );
[ Wh0 ] = Dirichlet_on_struct( Wh0 );

toc

fprintf('+--------------------------------------------------------------+\n');
fprintf('|               Start Riemannian Trust Regions...              |\n');
fprintf('+--------------------------------------------------------------+\n');

time_vect = zeros(number_of_runs,1);

for irun=1:number_of_runs
    fprintf('+--------------------------------------------------------------+\n');
    fprintf('| Run number: %d\n', irun );
    fprintf('+--------------------------------------------------------------+\n');
    % Solve the problem with Riemannian Trust Regions
    tic
    [ Wh, ~, info_tr, ~ ] = trustregions( problem, Wh0, options_tr );
    time_vect(irun) = toc;
    tr_gradnorm = [info_tr.gradnorm];
    num_inner_array = [info_tr.numinner];
    ntot_iter = length(tr_gradnorm) - 1;  % do not count the "0th" iteration
    fprintf('Total time: %.2f s.\n', time_vect(irun) );
    fprintf('Last grad. norm: %.4e.\n', tr_gradnorm(end) );
    fprintf('Residual: %.4e (norm of the Euclidean gradient).\n', norm_struct( problem.egrad( Wh ) ) );
    % CFR norm_Gh = norm_struct( getTruncationStruct( problem.egrad( Wh ), pars.K ) )
end

fprintf('****************************************************************\n');
fprintf('Number of outer iterations: %.2d\n', ntot_iter );
fprintf('Sum number of inner iters:  %.2d\n', sum(num_inner_array) );
fprintf('Max number of inner iters:  %.2d\n', max(num_inner_array) );

fprintf('****************************************************************');
Times = time_vect'
fprintf('  Number of runs: %d. Average time: %.2f s.\n', number_of_runs, mean(time_vect) );
fprintf('****************************************************************\n');

% diary off

% Plot
% Plot_PrecRTR( tr_gradnorm, pars );

% &&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%--------------------------------------------------------------------------
% SAVE DATA TO MAT-FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save data                          |\n');
fprintf('+--------------------------------------------------------------+\n');

fileName_mfile = [ 'mdatafiles/', pars.problem_type, '_Prec', num2str(pars.precon), ...
    '_Grad_', num2str(pars.lev_finest), '_K', num2str(pars.K), '.mat'];

save( fileName_mfile, 'tr_gradnorm', 'ntot_iter', 'Wh');
fprintf('Saved data to file %s.\n', fileName_mfile);
%--------------------------------------------------------------------------
