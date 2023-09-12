%==========================================================================
% Driver for plotting the time history of
%    \| W - W_{\mathrm{ref}} \|_{L^{2}(\Omega)}

% Created:     2023.04.10
% Last change: 2023.04.18

%   Apr 10, 2023:
%       Created by copying the original file.
%==========================================================================

close all; clear; clc;

options_plot;

% addpath(genpath('export_fig-master'));

%--------------------------------------------------------------------------
% gray3 = [ 0.255, 0.255, 0.255 ];
% super_light_gray = [ 0.99, 0.99, 0.99 ];

% % Added on 2023.03.06:
% set( 0, 'DefaultAxesXColor', super_light_gray );
% set( 0, 'DefaultAxesYColor', super_light_gray );
% set( 0, 'DefaultAxesZColor', super_light_gray );
% set( 0, 'DefaultTextColor', super_light_gray );
%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Spatial discretization:
% Interval boundaries:
wL = 0;
wR = 40;

% Spatial discretization:
Nx = 100;   % number of spatial grid points

x = linspace( wL, wR, Nx)';

% Mesh size
hx = (wR-wL)/(Nx-1);
%--------------------------------------------------------------------------
% Time evolution:
t0 = 0;   % initial time of the simulations
T = 1;     % final time of the simulations
%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------
% Reference solution:
load("../reference_solutions/FKPP_W_ref_hist_Nx_100_K_dt0.01.mat")
% load("../reference_solutions/FKPP_ERK4_W_ref_hist_Nx_100_K_dt0.0001.mat")
%--------------------------------------------------------------------------

% Load the saved data of the LRIE solution:
% load("../mdatafiles/FKPP_Prec1_W_hist_lra_Nx_100_K15_dt0.1.mat")
% load("../mdatafiles/FKPP_Prec1_W_hist_lra_Nx_100_K15_dt0.01.mat")

load("../reference_solutions/FKPP_Euclidean_W_final_100x100_dt0.01.mat")

% W_lra_mat = W_hist_lra(1).U * W_hist_lra(1).S * W_hist_lra(1).V';

% Colors in the yellow-orange palette
rnd_color_yellow_orange = [ ones(Nx,1), rand(Nx,1), 0.259*ones(Nx,1) ];

% [ U, S, V ] = svd( W_to_save_hist(:,:,1) );

% boolean_vect = diag(S) > tol_rank;

% Define the new rank:
% pars.K = nnz(boolean_vect);

% low_rank = 15;

% % Truncate
% W_ref_0.U = U(:,1:low_rank);
% W_ref_0.S = S(1:low_rank,1:low_rank);
% W_ref_0.V = V(:,1:low_rank);

% W_ref_0_truncated = W_ref_0.U*W_ref_0.S*W_ref_0.V';
% This should be the same as W_hist_lra(1).U * W_hist_lra(1).S * W_hist_lra(1).V'

% % some checks...
% norm( W_ref_0_truncated - W_to_save_hist(:,:,1) )
% norm( W_to_save_hist(:,:,1) - W_lra_mat )

norm( W_hist_Euclidean_TR(:,:,1) - W_ref_hist(:,:,1), "fro" ) * sqrt(hx)

figure(1)
% norm(W_ref_hist(:,:,end) - W_lra_mat, 'fro') / norm(W_lra_mat, "fro" )

% Plot numerical solution at final time:
plot( x, W_hist_Euclidean_TR(:,7,end), "LineWidth", 1 )
hold on
plot( x, W_ref_hist(:,7,end), "LineWidth", 1 )

% Random colors in the yellow-orange palette:
colororder( rnd_color_yellow_orange )
xlabel('$ x $')
ylabel('$ w $')
xlim( [ wL 25 ] )
% ylim( [ -0.02 1.02 ] )
title(['Low-rank+IE, time = ', num2str(T)])
drawnow;


% figure(2)
% % Plot a single trajectory at the final time:
% i = 3;
% plot( x, W_ref_hist(:,i,end), '--', 'LineWidth', 2, 'Color', blue )
% hold on
% plot( x, W_lra_mat(:,i), 'LineWidth', 2, 'Color', red )
% grid on
% xlabel('$ x $')
% ylabel('$ w $')
% legend( 'reference solution for trajectory $ w^{(500)} $', ...
%     'LRIE solution for trajectory $ w^{(500)} $' )
% xlim( [ uL uR ] )
% ylim( [ -0.03 1.01 ] )

% norm(U_t_hist(:,i,end)-W_hist_lra_sel(:,i))*sqrt(hx)

% return
% 
for i=1:size(W_hist_Euclidean_TR,3)
%     W_hist_lra_sel = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt(i) = norm( W_hist_Euclidean_TR(:,:,i) - W_ref_hist(:,:,i), "fro" ) * sqrt(hx);
%     L2_norm_dt_first_trajectory(i) = norm(U_t_hist(:,1,end)-W_hist_lra_sel(:,1)) * sqrt(hx);
end


%--------------------------------------------------------------------------
% Plot
figure(3)
semilogy( t_hist, L2_norm_dt, "Color", red, "LineStyle", "--", "LineWidth", 2 )
grid on
%--------------------------------------------------------------------------
% MS, 2023.03.06:
% Change the color of the grid:
% set( gca,'GridColor', super_light_gray )

% Added on 2023.03.06:
% Questo è il background "esterno" alla figura
% set( gca,'Color', gray3 )
%--------------------------------------------------------------------------
xlabel('$t$')
ylabel('$ \| w - w_{\mathrm{ref}} \|_{L^{2}(\Omega)} $')
legend( '$ h = 0.125 $' )
xlim( [ t_hist(1), t_hist(end) ] )
ylim( [ 9e-10 9e1 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                         Save figure 1                        |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = '../plots/FKPP_L2_norm_W_lra_vs_W_ref_rank_adaptive';

set(gcf, 'Color', 'w');
export_fig(fileName_plot, '-pdf', '-opengl','-r600');

% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.pdf.\n', fileName_plot);
