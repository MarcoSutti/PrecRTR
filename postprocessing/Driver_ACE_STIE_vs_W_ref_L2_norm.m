%==========================================================================
% Driver for plotting the time evolution of 
%    \| W - W_{\mathrm{ref}} \|_{L^{2}(\Omega)}
% for RV2022 Step-Truncation Implicit Euler.

% Created:     2023.02.17
% Last change: 2023.02.22

%   Feb 17, 2023:
%       Created by copying the original file.
%==========================================================================

close all; clear; clc;

options_plot;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Spatial discretization:
Nx = 256;
Lx = 2*pi;
hx = Lx/Nx;
%--------------------------------------------------------------------------
% Time evolution:
% t0 = 0.5;   % initial time of the simulations
T = 15;     % final time of the simulations
%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------
% Reference solution:
load("ACE_ref_256x256_T25_dt0.0001.mat")

% Cut reference solution at between t0 and T:
t_hist_ref = t_hist(1:151);
W_hist_ref = W_hist(:,:,1:151);
%--------------------------------------------------------------------------
% Load the saved data:
load("ACE_STIE_256x256_T15_dt0.0005_from_t0_0.mat")

for i=1:length(t_STIE_hist)
    L2_norm_dt1(i) = norm( W_hist_ref(:,:,i) - W_STIE_hist(:,:,i), 'fro' ) * hx;
end


%--------------------------------------------------------------------------
% Plot
figure(1)
semilogy( t_STIE_hist, L2_norm_dt1, "Color", blue, "LineStyle", "-", "LineWidth", 2 )
grid on
xlabel('$t$')
ylabel('$ \| w - w_{\mathrm{ref}} \|_{L^{2}(\Omega)} $')
legend( '$ h = 10^{-3} $' )
xlim( [ 0 15 ] )
ylim( [ 9e-10 9e1 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                         Save figure 1                        |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = 'plots/ACE_L2_norm_W_STIE_vs_W_ref_256x256';

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName_plot);
%--------------------------------------------------------------------------
%%
% % Plot the "stability barrier"
% figure(2)
% dt_array = [ 0.05, 0.1, 0.2, 0.5, 1.0 ];
% L2_norm_final_T = [ L2_norm_dt2(end), L2_norm_dt1(end), ...
%     L2_norm_dt5_T15, L2_norm_dt3(end), L2_norm_dt4_T15 ];
% 
% loglog( dt_array, L2_norm_final_T, 'o-', 'Color', red, ...
%     'LineWidth', 2, ...
%     'MarkerFaceColor', red, 'MarkerSize', 8 )
% hold on
% loglog( dt_array, 5e-5 * dt_array.^1, '--', 'Color', gray2 )
% grid on
% legend( 'error', '$ \mathcal{O}(h) $', 'FontSize', 14, 'Location', 'SE' )
% xlabel('$ h $')
% ylabel('$ \| w_{T} - w_{\mathrm{ref},T} \|_{L^{2}(\Omega)} $')
% ylim( [ 1e-6 5e-4 ] )
% 
% %--------------------------------------------------------------------------
% % SAVE IMAGE TO EPS FILE
% %--------------------------------------------------------------------------
% fprintf('+--------------------------------------------------------------+\n');
% fprintf('|                         Save figure 2                        |\n');
% fprintf('+--------------------------------------------------------------+\n');
% fileName_plot = [ 'plots/ACE_stability_barrier_t0_', num2str(t0) ];
% 
% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
% fprintf('Saved graph to file %s.eps.\n', fileName_plot);
% %--------------------------------------------------------------------------

% Plot
figure(2)
semilogy( t_STIE_hist, rank_STIE_hist, "Color", red, "LineStyle", "-", "LineWidth", 2 )
grid on
xlabel('$t$')
ylabel('$ \mathrm{rank}(W_{\mathrm{ref}}) $')
% xlim( [ 0, 3 ] )
% ylim( [ min(rank_STIE_hist)-2  max(rank_STIE_hist)+20 ] )
ylim( [ 1, 100 ])