%==========================================================================
% Driver for 1D stochastic Fisher-KPP PDE with homogeneous Neumann BCs.
% This script compares the solution given by the IMEX-CNLF method to the
% solution provided by the Euclidean optimization on an appropriately built
% cost function.
% Multiple realizations.

% Created:     2023.04.27
% Last change: 2023.05.08

%   May 8, 2023:
%       Cleanup of the code, and added save of CNLF reference solution.
%   Apr 27, 2023:
%       Created.
%==========================================================================

close all; clear; clc;

% Fixed rng seed:
rng(1);

addpath( genpath('plots') );
addpath( genpath('utilities') );
addpath( genpath('reference_solutions') );

options_plot;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Interval boundaries:
wL = 0;
wR = 40;

% Spatial discretization:
Nx = 1000;   % number of spatial grid points
x = linspace( wL, wR, Nx)';

% Total number of realizations:
Nr = Nx;
%--------------------------------------------------------------------------
% Time discretization
T = 10;   % final time
Nt = 401; % points in time
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Parameters for Euclidean trustregions (Manopt)
options_EuTR.maxiter = 100;
options_EuTR.minstepsize = 1e-12;
options_EuTR.tolgradnorm = 1e-10;
% options_tr.Delta_bar = 1;     % MS, 15.01.2021: It is better to comment
% %                               out this line and let manopt take care of
% %                               the Delta_bar automatically.
options_EuTR.verbosity = 0;
%--------------------------------------------------------------------------
% Parameters for Riemannian trustregions (Manopt)
options_RTR.maxiter = 100;
options_RTR.minstepsize = 1e-12;
options_RTR.tolgradnorm = 1e-8;
options_RTR.Delta_bar = 1;
options_RTR.verbosity = 0;

% Preconditiong: 0 for No, 1 for Yes.
pars.precon = 1;

% Manifold parameters
% Retraction type: 'metric', 'metric_explicit_inverse', 'ortho_ao', 'ortho_zhang'.
pars.retr_type = 'metric';

% Tolerance for rank truncation:
pars.tol_rank = 1e-8;

%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Spatial discretization
hx = (wR-wL)/(Nx-1);

%--------------------------------------------------------------------------
% Save in the pars structure
pars.hx = hx;
pars.Nx = Nx;
pars.Nr = Nr;
pars.Nt = Nt;
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                          Precompute A                        |\n');
fprintf('+--------------------------------------------------------------+\n');
% Create matrix A in sparse format
% The discretized Laplacian with homogeneous Neumann boundary conditions:
Ah = get_Ah_Neumann_BCs( Nx, hx );
pars.A = Ah;
%--------------------------------------------------------------------------
% Time step
h = T/(Nt-1);

% Save in the pars struct:
pars.dt = h;

%--------------------------------------------------------------------------
pars.Mplus = speye(Nx) + pars.dt * pars.A;
pars.Mminus = speye(Nx) - pars.dt * pars.A;

pars.MmtMm = pars.Mminus'*pars.Mminus;
% pars.MmtMm = speye(Nx) - pars.dt*pars.A' - pars.dt*pars.A + pars.dt^2 * pars.A'*pars.A;
pars.MptMm = pars.Mplus'*pars.Mminus;
%--------------------------------------------------------------------------

t_hist = zeros(1, Nt);

fprintf('+--------------------------------------------------------------+\n');
fprintf('|                  Define initial conditions...                |\n');
fprintf('+--------------------------------------------------------------+\n');
%--------------------------------------------------------------------------
% 1. Define initial conditions and use ERK4 for the first time step.
%--------------------------------------------------------------------------
[ W0, Romega ] = get_FKPP_IC( Nx, Nr, x );

time_iter = 1;
current_t = (time_iter-1) * h;
fprintf( "IC: Time: %.4f \n", current_t );

pars.Romega = Romega;

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%--------------------------------------------------------------------------
% SAVE INITIAL CONDITIONS TO MAT-FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                   Save initial conditions                    |\n');
fprintf('+--------------------------------------------------------------+\n');

% W0 = W_CNLF_hist(:,:,1);

fileName_mfile = [ 'reference_solutions/FKPP_ICs_Nx', num2str(Nx), ...
    '_Nr', num2str(Nr) ];

% Save W0 and Romega:
save( fileName_mfile, 'W0', 'Romega' );
fprintf('Saved data to file %s.mat.\n', fileName_mfile);
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 1.2. Use ERK4 for the first time step.
%      NB: For all the realizations at once. See my notes of 2023.04.27.
%--------------------------------------------------------------------------
time_iter = 2;
current_t = (time_iter-1) * h;
fprintf( "Time: %.4f \n", current_t );

% Define the function on the right-hand side:
fun_RHS = @(W) -pars.A * W + W.*(1-W)*pars.Romega;

Wn_minus_1 = one_step_ERK4( W0, h, fun_RHS );

% Stores the first approximate solution computed:
t_hist(time_iter) = current_t;
fprintf( "ERK4: Time: %.4f \n", current_t );

%--------------------------------------------------------------------------
% 2. Now use the IMEX-CNLF method or the Euclidean TR method for all the
%    other time steps.
%    We treat all the realizations at once. See my notes of 2023.04.27.
[ t_EuTR_hist, W_EuTR_hist, t_hist_stride, W_CNLF_hist_stride, rank_W_CNLF_hist_stride ] = integrate_IMEX_CNLF( W0, Wn_minus_1, pars, options_EuTR );

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%--------------------------------------------------------------------------
% SAVE SOLUTION AT EVERY STRIDE T TO MAT-FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|       Save IMEX-CNLF numerical solution at every stride      |\n');
fprintf('+--------------------------------------------------------------+\n');

fileName_mfile = [ 'reference_solutions/FKPP_CNLF_W_hist_Nx', num2str(Nx), ...
    '_Nr', num2str(Nr), '_T', num2str(T), '_Nt', num2str(Nt) ];

save( fileName_mfile, 't_hist_stride', 'W_CNLF_hist_stride', 'rank_W_CNLF_hist_stride' );
fprintf('Saved data to file %s.mat.\n', fileName_mfile);
%--------------------------------------------------------------------------
%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&


% Here you can change the time parameters!!!!!!!!!!!!!!!!!!!!!!!!!!
Nt = pars.Nt; % points in time

% Time step
h = T/(Nt-1);

% Save in the pars struct:
pars.dt = h;

%--------------------------------------------------------------------------
pars.Mplus = speye(Nx) + pars.dt * pars.A;
pars.Mminus = speye(Nx) - pars.dt * pars.A;

pars.MmtMm = pars.Mminus'*pars.Mminus;
pars.MptMm = pars.Mplus'*pars.Mminus;
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 2.2. Use ERK4 for the first time step.
%      NB: For all the realizations at once. See my notes of 2023.04.27.
%--------------------------------------------------------------------------
time_iter = 2;
current_t = (time_iter-1) * h;
%     fprintf( "Time: %.4f \n", current_t );

% !!! NB: r_omega changes for each realization !!!
fun_RHS = @(W) -pars.A * W + W.*(1-W)*pars.Romega;

Wn_minus_1 = one_step_ERK4( W0, h, fun_RHS );

% Stores the first approximate solution computed:
t_hist(time_iter) = current_t;
fprintf( "ERK4: Time: %.4f \n", current_t );

%--------------------------------------------------------------------------
% 2.3 Riemannian Trust-Region Method
%--------------------------------------------------------------------------
[ Wn_RTR_hist, rank_W_RTR_hist ] = integrate_LR_CNLF( W0, Wn_minus_1, pars, options_RTR );

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%--------------------------------------------------------------------------
% SAVE LR-CNLF SOLUTION TO MAT-FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                Save LR-CNLF numerical solution               |\n');
fprintf('+--------------------------------------------------------------+\n');

fileName_mfile = [ 'reference_solutions/FKPP_RTR_W_hist_Nx', num2str(Nx), ...
    '_Nr', num2str(Nr), '_T', num2str(T), '_Nt', num2str(Nt), ...
    '_final_rank', num2str(rank_W_RTR_hist(end)) ];

save( fileName_mfile, 'Wn_RTR_hist', 'rank_W_RTR_hist' );
fprintf('Saved data to file %s.mat.\n', fileName_mfile);
%--------------------------------------------------------------------------



%%
clc
% W_CNLF_hist_stride

stride = floor(size(Wn_RTR_hist,2)/100);

norm_diff_W = zeros(1,size(W_CNLF_hist_stride,3));

k = 0;
for i=1:stride:size(Wn_RTR_hist,2)
    k = k + 1;
    norm_diff_W(k) = norm( Wn_RTR_hist(i).U * Wn_RTR_hist(i).S * Wn_RTR_hist(i).V' - W_CNLF_hist_stride(:,:,k), 'fro' ) * sqrt(hx);
end


figure(1)
close all

semilogy( t_EuTR_hist(1:stride:size(Wn_RTR_hist,2)), norm_diff_W, "Color", red, "LineStyle", "--", "LineWidth", 3 )
grid on