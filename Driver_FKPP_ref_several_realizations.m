%==========================================================================
% Driver for 1D stochastic Fisher-KPP PDE with homogeneous Neumann BCs.
% From the paper of Charous and Lermusiaux, 2022.
% NB: This script runs the simulation. To plot the generated data, please
%     use the script 'Driver_plot_FKPP_ref_several_realizations'.
% This is a reproduction of the results in the preprint:
% [1] Stable rank-adaptive dynamically orthogonal Runge-Kutta schemes,
%     Charous, A. and Lermusiaux, P., available on arXiv, 2022.

% Created:     2023.03.20
% Last change: 2023.04.25

%   Apr 25, 2023:
%       Modified the way the reference solution is stored (ns).
%   Mar 22, 2023:
%       Created by copying Driver_Fisher_KPP.
%==========================================================================

close all; clear; clc;

% Fixed rng seed:
rng(1);

addpath( genpath('plots') );
addpath( genpath('utilities') );
addpath( genpath('reference_solutions') );

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Interval boundaries:
wL = 0;
wR = 40;

% Spatial discretization:
Nx = 1000;   % number of spatial grid points
x = linspace( wL, wR, Nx)';

% Time discretization
T = 12.5;   % final time
Nt = 10001; % points in time

% Total number of realizations:
tot_num_realizations = 1000;

% Save data every ns iterations
ns = 100;
%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Spatial discretization
hx = (wR-wL)/(Nx-1);
%--------------------------------------------------------------------------
% The discretized minus Laplacian (without the prefactor 1/h^2) in sparse
% format with Neumann BCs:
A = get_Ah_Neumann_BCs( Nx );
%--------------------------------------------------------------------------
% Time step
h = T/(Nt-1);
%--------------------------------------------------------------------------
% Ratio of dt over dx squared:
r = h/hx^2;   % the best would be to have this quantity neither too
              % big nor too small, so that the matrix M is well conditioned
%--------------------------------------------------------------------------

M1 = speye(Nx) + r * A;
M2 = speye(Nx) - r * A;

fprintf( "Cond number of M1: %.f \n", condest(M1) );
fprintf( "Cond number of M2: %.f \n", condest(M2) );
%--------------------------------------------------------------------------

% Initializations:
r_omega_array = zeros(1, tot_num_realizations);
t_hist = zeros(1, Nt);
% Dummy variables used in the for loops:
W_minus_1 = zeros(Nx, tot_num_realizations);
W_minus_2 = zeros(Nx, tot_num_realizations);

Nt_ns = (Nt-1)/ns + 1;
W_ns_hist = zeros(Nx, tot_num_realizations, Nt_ns);
t_ns_hist = zeros(1, Nt_ns);

for idx_realization=1:tot_num_realizations
    fprintf( "Realization: %.d \n", idx_realization );
    %----------------------------------------------------------------------
    % Define initial condition.
    % The coefficients a, b, r follow a uniform law.
    a_omega = unifrnd(1/5, 2/5);
    b_omega = unifrnd(1/10, 11/10);
    r_omega = unifrnd(1/4, 1/2);   % reaction rate

    % Initial condition:
    w0 = a_omega * exp(-b_omega * x.^2);

    % Store u0 in U_0 for later saving.
    time_iter = 1;
    current_t = (time_iter-1) * h;
    t_hist(time_iter) = current_t;
    W_minus_2(:,idx_realization) = w0;
    
    % Save initial condition:
    idx_hist = 1;
    W_ns_hist(:,idx_realization,idx_hist) = w0;
    t_ns_hist(idx_hist) = current_t;

    % Store r_omega in r_omega_array for later saving
    r_omega_array(idx_realization) = r_omega;
    %----------------------------------------------------------------------

    r2 = 2*h*r_omega;

    %----------------------------------------------------------------------
    % Time evolution.
    % Use ERK4 for the first step.
    time_iter = 2;
    current_t = (time_iter-1) * h;

    % !!! NB: r_omega changes for each realization !!!
    fun_rhs = @(u) (1/hx^2) * A * u + r_omega * u.*(1-u);

    k1 = fun_rhs( w0 );
    k2 = fun_rhs( w0 + (h/2)*k1 );
    k3 = fun_rhs( w0 + (h/2)*k2 );
    k4 = fun_rhs( w0 + h*k3 );
    wn_minus_1 = w0 + (h/6)*( k1 + 2*k2 + 2*k3 + k4 );

    % Stores the first approximate solution computed:
    W_minus_1(:,idx_realization) = wn_minus_1;
    t_hist(time_iter) = current_t;

    % Now use the IMEX-CNLF method for all the other time steps.
    for time_iter=3:Nt
        wn_minus_1 = W_minus_1(:,idx_realization);
        wn_minus_2 = W_minus_2(:,idx_realization);

        b = M2*wn_minus_2 + r2*wn_minus_1.*(1-wn_minus_1);

        wn = M1\b;

        % Store the new approximate solution:
        current_t = (time_iter-1) * h;
        t_hist(time_iter) = current_t;
        W_minus_2(:,idx_realization) = W_minus_1(:,idx_realization);
        W_minus_1(:,idx_realization) = wn;

        % fprintf('%d  \t %.3f \n', time_iter, current_t );
        if mod(time_iter-1,ns)==0
            idx_hist = idx_hist + 1;
            W_ns_hist(:,idx_realization,idx_hist) = wn;
            t_ns_hist(idx_hist) = current_t;
        end

    end

end


%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%--------------------------------------------------------------------------
% SAVE DATA TO MAT-FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save data                          |\n');
fprintf('+--------------------------------------------------------------+\n');

fileName_mfile = [ 'reference_solutions/FKPP_W_ref_hist_Nx_', num2str(Nx), ...
    '_K_dt', num2str(h), '.mat'];

save( fileName_mfile, 'W_ns_hist', 't_ns_hist', 'r_omega_array' )
fprintf('Saved data to file %s.\n', fileName_mfile);
