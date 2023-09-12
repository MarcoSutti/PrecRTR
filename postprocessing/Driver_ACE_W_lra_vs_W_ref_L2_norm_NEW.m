%==========================================================================
% Driver for plotting the time evolution of 
%    \| W - W_{\mathrm{ref}} \|_{L^{2}(\Omega)}

% Created:     2023.02.13
% Last change: 2023.03.06

%   Mar 6, 2023:
%       Added lines of codes to save a transparent png version of the plot,
%       in order to generate the plot for the talk in the NCTS Spring Day.
%   Feb 14, 2023:
%       Replaced arithmetic mean with geometric mean.
%   Feb 13, 2023:
%       Created by copying the original file.
%==========================================================================

close all; clear; clc;

options_plot;

addpath(genpath('export_fig-master'));

%--------------------------------------------------------------------------
gray3 = [ 0.255, 0.255, 0.255 ];
super_light_gray = [ 0.99, 0.99, 0.99 ];

% Added on 2023.03.06:
set( 0, 'DefaultAxesXColor', super_light_gray );
set( 0, 'DefaultAxesYColor', super_light_gray );
set( 0, 'DefaultAxesZColor', super_light_gray );
set( 0, 'DefaultTextColor', super_light_gray );
%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Spatial discretization:
Nx = 256;
Lx = 2*pi;
hx = Lx/Nx;
%--------------------------------------------------------------------------
% Time evolution:
t0 = 0.5;   % initial time of the simulations
T = 15;     % final time of the simulations
%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------
% Reference solution:
load("ACE_ref_256x256_T25_dt0.0001.mat")

% Cut reference solution at between t0 and T:
t_hist_ref = t_hist(6:151);
W_hist_ref = W_hist(:,:,6:151);
%--------------------------------------------------------------------------
% dt = 0.1
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K11_dt0.1.mat")

t_hist_lra_1 = t_hist_lra(1:146);

% offset = length(t_hist_ref)-length(t_hist_lra_1);

for i=1:length(t_hist_ref)
    W_hist_lra_1 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt1(i) = norm( W_hist_ref(:,:,i) - W_hist_lra_1, 'fro' ) * hx;
end

% %--------------------------------------------------------------------------
% % dt = 0.01
% % Load the saved data:
% load("ACE_Prec1_W_hist_lra_256x256_K13_dt0.01.mat")
% 
% select_stride_2 = 1:10:length(t_hist_lra);
% t_hist_lra_2 = t_hist_lra(select_stride_2);
% 
% k = 0;
% for i=1:10:size(W_hist_lra,2)
%     k = k + 1;
%     W_hist_lra_2(:,:,k) = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
% end
% 
% offset = length(t_hist_ref)-length(t_hist_lra_2);
% 
% for i=1:length(t_hist_lra_2)
%     L2_norm_dt2(i) = norm( W_hist(:,:,i + offset) - W_hist_lra_2(:,:,i), 'fro' ) * hx;
% end

%--------------------------------------------------------------------------
% dt = 0.05
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K11_dt0.05.mat")

select_stride_2 = 1:2:291;
t_hist_lra_2 = t_hist_lra(select_stride_2);

k = 0;
for i=select_stride_2
    k = k + 1;
    W_hist_lra_2(:,:,k) = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
end

for i=1:length(t_hist_ref)
    L2_norm_dt2(i) = norm( W_hist_ref(:,:,i) - W_hist_lra_2(:,:,i), 'fro' ) * hx;
end

%--------------------------------------------------------------------------
% dt = 0.20
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K17_dt0.2.mat")

% select_stride_5 = 1:2:length(t_hist_ref);
% t_hist_5 = t_hist_ref(select_stride_5);
% W_hist_ref_5 = W_hist_ref(:,:,select_stride_5);
W_hist_ref_5 = W_hist(:,:,6:2:153);

t_hist_lra_5 = t_hist_lra;

for i=1:length(t_hist_lra_5)
    W_hist_lra_5 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt5(i) = norm( W_hist_ref_5(:,:,i) - W_hist_lra_5, 'fro' ) * hx;
end

%--------------------------------------------------------------------------
% dt = 0.50
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K12_dt0.5.mat")

select_stride_3 = 1:5:length(t_hist_ref);
% t_hist_3 = t_hist_ref(select_stride_3);
W_hist_ref_3 = W_hist_ref(:,:,select_stride_3);
W_hist_ref_4 = W_hist(:,:,6:10:161);

t_hist_lra_3 = t_hist_lra(1:end-1);

for i=1:length(t_hist_lra_3)
    W_hist_lra_3 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt3(i) = norm( W_hist_ref_3(:,:,i) - W_hist_lra_3, 'fro' ) * hx;
end


%--------------------------------------------------------------------------
% dt = 1.00
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K14_dt1.mat")

% select_stride_4 = 1:10:length(t_hist_ref);
% t_hist_4 = t_hist_ref(select_stride_4);
W_hist_ref_4 = W_hist(:,:,6:10:161);

t_hist_lra_4 = t_hist_lra;

for i=1:length(t_hist_lra_4)
    W_hist_lra_4 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt4(i) = norm( W_hist_ref_4(:,:,i) - W_hist_lra_4, 'fro' ) * hx;
end

% Linear interp. to get the L2-norm at T=15 for the case dt = 1.00
% L2_norm_dt4_T15 = 0.5 * ( L2_norm_dt4(end-1) + L2_norm_dt4(end) );
% Take the geometric mean instead of the arithmetic mean because we are
% using the semilogarithmic scale:
L2_norm_dt4_T15 = sqrt( L2_norm_dt4(end-1) * L2_norm_dt4(end) );

% Geometric mean to get the L2-norm at T=15 for the case dt = 0.20
L2_norm_dt5_T15 = sqrt( L2_norm_dt5(end-1) * L2_norm_dt5(end) );

%--------------------------------------------------------------------------
% Plot
figure(1)
semilogy( t_hist_lra_4, L2_norm_dt4, "Color", red, "LineStyle", "--", "LineWidth", 2 )
hold on
semilogy( t_hist_lra_3, L2_norm_dt3, "Color", blue, "LineStyle", "-.", "LineWidth", 2 )
semilogy( t_hist_lra_5, L2_norm_dt5, "Color", orange, "LineStyle", "-", "LineWidth", 2 )
semilogy( t_hist_lra_1, L2_norm_dt1, "Color", green, "LineStyle", "-", "LineWidth", 2 )
semilogy( t_hist_lra_2, L2_norm_dt2, "Color", yellow, "LineStyle", "--", "LineWidth", 2 )
grid on
%--------------------------------------------------------------------------
% MS, 2023.03.06:
% Change the color of the grid:
set( gca,'GridColor', super_light_gray )

% Added on 2023.03.06:
% Questo è il background "esterno" alla figura
set( gca,'Color', gray3 )
%--------------------------------------------------------------------------
xlabel('$t$')
ylabel('$ \| w - w_{\mathrm{ref}} \|_{L^{2}(\Omega)} $')
legend( '$ h = 1.0 $', '$ h = 0.5 $', '$ h = 0.20 $', '$ h = 0.1 $', '$ h = 0.05 $' )
% semilogy( 15, L2_norm_dt4_T15, 'o-', 'Color', red, ...
%     'LineWidth', 2, ...
%     'MarkerFaceColor', red, 'MarkerSize', 8 )
% semilogy( 15, L2_norm_dt5_T15, 'x-', 'Color', 'k', ...
%     'LineWidth', 2, ...
%     'MarkerFaceColor', red, 'MarkerSize', 8 )
xlim( [ 0 15 ] )
ylim( [ 9e-10 9e1 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                         Save figure 1                        |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = [ 'plots/ACE_L2_norm_W_lra_vs_W_ref_256x256_rank_adaptive_t0_', num2str(t0) ];

% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
% fprintf('Saved graph to file %s.eps.\n', fileName_plot);
% %--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Save plot to png file
export_fig(fileName_plot, '-dpng', '-transparent', '-r300'); % -r300 is the PPI value, default resolution is low
fprintf('Saved graph to file %s.png.\n', fileName_plot);
%--------------------------------------------------------------------------
%%
% Plot the "stability barrier"
figure(2)
dt_array = [ 0.05, 0.1, 0.2, 0.5, 1.0 ];
L2_norm_final_T = [ L2_norm_dt2(end), L2_norm_dt1(end), ...
    L2_norm_dt5_T15, L2_norm_dt3(end), L2_norm_dt4_T15 ];

loglog( dt_array, L2_norm_final_T, 'o-', 'Color', red, ...
    'LineWidth', 2, ...
    'MarkerFaceColor', red, 'MarkerSize', 8 )
hold on
loglog( dt_array, 5e-5 * dt_array.^1, '--', 'Color', gray2 )
grid on
%--------------------------------------------------------------------------
% MS, 2023.03.06:
% Change the color of the grid:
set( gca,'GridColor', super_light_gray )

% Added on 2023.03.06:
% Questo è il background "esterno" alla figura
set( gca,'Color', gray3 )
%--------------------------------------------------------------------------
legend( 'error', '$ \mathcal{O}(h) $', 'FontSize', 14, 'Location', 'SE' )
xlabel('$ h $')
ylabel('$ \| w_{T} - w_{\mathrm{ref},T} \|_{L^{2}(\Omega)} $')
ylim( [ 1e-6 5e-4 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                         Save figure 2                        |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = [ 'plots/ACE_stability_barrier_t0_', num2str(t0) ];

% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
% fprintf('Saved graph to file %s.eps.\n', fileName_plot);
% %--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% Save plot to png file
export_fig(fileName_plot, '-dpng', '-transparent', '-r300'); % -r300 is the PPI value, default resolution is low
fprintf('Saved graph to file %s.png.\n', fileName_plot);
%--------------------------------------------------------------------------