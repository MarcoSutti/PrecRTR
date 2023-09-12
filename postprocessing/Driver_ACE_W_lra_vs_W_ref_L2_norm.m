%==========================================================================
% Driver for plotting the time evolution of 
%    \| W - W_{\mathrm{ref}} \|_{L^{2}(\Omega)}

% Created:     2023.01.10
% Last change: 2023.02.04

%   Feb 4, 2023:
%       Tried to compare the errors at t = 8 instead of t = T, it seems
%       that the error decays linearly with dt.
%       Added solution with dt = 0.001, and dt = 0.005. 
%   Jan 29, 2023:
%       Added solution with dt = 0.5.
%   Jan 11, 2023:
%       Added solution with dt = 1e-2.
%   Jan 10, 2023:
%       Added solution with the CFD discretization of the Laplacian.
%       Created.
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
% End of data
%--------------------------------------------------------------------------

% dt = 1e-3
% Load the saved data:
load("ACE_ref_T25_dt0.0001.mat")
load("ACE_Prec1_W_hist_lra_256x256_K10_dt0.001bis.mat")


select_stride_1 = 1:10:length(t_hist_lra);
t_hist_lra_1 = t_hist_lra(select_stride_1);

offset = length(t_hist)-length(t_hist_lra_1);

for i=1:length(t_hist_lra_1)
    W_hist_lra_1 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt1(i) = norm( W_hist(:,:,i + offset) - W_hist_lra_1, 'fro' ) * hx;
end

%--------------------------------------------------------------------------
% dt = 5e-3
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K10_dt0.005bis.mat")


select_stride_2 = 1:2:length(t_hist_lra);
t_hist_lra_2 = t_hist_lra(select_stride_2);

offset = length(t_hist)-length(t_hist_lra_2);

for i=1:length(t_hist_lra_2)
    W_hist_lra_2 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt2(i) = norm( W_hist(:,:,i + offset) - W_hist_lra_2, 'fro' ) * hx;
end


%--------------------------------------------------------------------------
% dt = 1e-2
% Load the saved data:
% load("ACE_ref_sol_Ah_T15_dt1e-2.mat")
load("ACE_Prec1_W_hist_lra_256x256_K11_dt0.01bis.mat")

t_hist_lra_3 = t_hist_lra;

offset = length(t_hist)-length(t_hist_lra_3);

for i=1:length(t_hist_lra_3)
    W_hist_lra_3 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt3(i) = norm( W_hist(:,:,i + offset) - W_hist_lra_3, 'fro' ) * hx;
end


%--------------------------------------------------------------------------
% dt = 5e-2
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K14_dt0.05bis.mat")

select_stride_4 = 1:5:length(t_hist);
t_hist_4 = t_hist(select_stride_4);
W_hist_4 = W_hist(:,:,select_stride_4);

t_hist_lra_4 = t_hist_lra;

offset = length(t_hist_4)-length(t_hist_lra);

for i=1:length(t_hist_lra)
    W_hist_lra_4 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt4(i) = norm( W_hist_4(:,:,i + offset) - W_hist_lra_4, 'fro' ) * hx;
end

%--------------------------------------------------------------------------
% dt = 0.1
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K12_dt0.1bis.mat")

select_stride_5 = 1:10:length(t_hist);
t_hist_5 = t_hist(select_stride_5);
W_hist_5 = W_hist(:,:,select_stride_5);

t_hist_lra_5 = t_hist_lra;

offset = length(t_hist_5)-length(t_hist_lra);

for i=1:length(t_hist_lra)
    W_hist_lra_5 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt5(i) = norm( W_hist_5(:,:,i + offset) - W_hist_lra_5, 'fro' ) * hx;
end

%--------------------------------------------------------------------------
% dt = 0.25
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K15_dt0.25bis.mat")

select_stride_6 = 1:25:length(t_hist);
t_hist_6 = t_hist(select_stride_6);
W_hist_6 = W_hist(:,:,select_stride_6);

t_hist_lra_6 = t_hist_lra;

offset = length(t_hist_6)-length(t_hist_lra);

for i=1:length(t_hist_lra)
    W_hist_lra_6 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt6(i) = norm( W_hist_6(:,:,i + offset) - W_hist_lra_6, 'fro' ) * hx;
end

%--------------------------------------------------------------------------
% dt = 0.50
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K11_dt0.5bis.mat")

select_stride_7 = 1:50:length(t_hist);
t_hist_7 = t_hist(select_stride_7);
W_hist_7 = W_hist(:,:,select_stride_7);

t_hist_lra_7 = t_hist_lra;

offset = length(t_hist_7)-length(t_hist_lra);

for i=1:length(t_hist_lra)
    W_hist_lra_7 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt7(i) = norm( W_hist_7(:,:,i + offset) - W_hist_lra_7, 'fro' ) * hx;
end


%--------------------------------------------------------------------------
% dt = 1.00
% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K12_dt1bis.mat")

select_stride_8 = 1:100:length(t_hist);
t_hist_8 = t_hist(select_stride_8);
W_hist_8 = W_hist(:,:,select_stride_8);

t_hist_lra_8 = t_hist_lra;

offset = length(t_hist_8)-length(t_hist_lra);

for i=1:length(t_hist_lra)
    W_hist_lra_8 = W_hist_lra(i).U * W_hist_lra(i).S * W_hist_lra(i).V';
    L2_norm_dt8(i) = norm( W_hist_8(:,:,i + offset) - W_hist_lra_8, 'fro' ) * hx;
end

%--------------------------------------------------------------------------
% Plot
figure(1)
semilogy( t_hist_lra_8, L2_norm_dt8, "Color", red, "LineStyle", ":", "LineWidth", 2 )
hold on
semilogy( t_hist_lra_7, L2_norm_dt7, "Color", red, "LineStyle", "--", "LineWidth", 2 )
semilogy( t_hist_lra_6, L2_norm_dt6, "Color", green, "LineStyle", "-.", "LineWidth", 2 )
semilogy( t_hist_lra_5, L2_norm_dt5, "Color", blue, "LineStyle", ":", "LineWidth", 2 )
semilogy( t_hist_lra_4, L2_norm_dt4, "Color", yellow, "LineStyle", "-.", "LineWidth", 2 )
semilogy( t_hist_lra_3, L2_norm_dt3, "Color", green, "LineStyle", "-", "LineWidth", 2 )
semilogy( t_hist_lra_2, L2_norm_dt2, "Color", blue, "LineStyle", "-.", "LineWidth", 2 )
semilogy( t_hist_lra_1, L2_norm_dt1, "Color", red, "LineStyle", ":", "LineWidth", 2 )
grid on
xlabel('$t$')
ylabel('$ \| w - w_{\mathrm{ref}} \|_{L^{2}(\Omega)} $')
legend( '$ \Delta t = 1.0 $', '$ \Delta t = 0.5 $', '$ \Delta t = 0.25 $', ...
    '$ \Delta t = 0.1 $', '$ \Delta t = 0.05 $', '$ \Delta t = 0.01 $', ...
    '$ \Delta t = 0.005 $', '$ \Delta t = 0.001 $' )
% xlim( [ 10, 15 ] )
ylim( [ 9e-10 9e1 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                         Save figure 1                        |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = 'plots/ACE_L2_norm_W_lra_vs_W_ref_256x256_rank_adaptive_from_t0.1';

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName_plot);
%--------------------------------------------------------------------------
%%
% Plot the "stability barrier"
figure(2)
% % dt_array = [ 0.01, 0.05, 0.1, 0.25, 0.5, 1.0 ];
% dt_array = [ 0.1, 0.5, 1.0 ];
% % L2_norm_final_T = [ L2_norm_dt1(end), L2_norm_dt4(end), ...
% %     L2_norm_dt2(end), L2_norm_dt5(end), L2_norm_dt3(end), L2_norm_dt6(end) ];
% L2_norm_final_T = [ L2_norm_dt2(end), L2_norm_dt3(end), L2_norm_dt6(end) ];
% 
% % At t = 10
% dt_array = [ 0.1, 0.25, 0.5, 1.0 ];
% L2_norm_final_T = [ L2_norm_dt2(96), L2_norm_dt5(39), L2_norm_dt3(20), L2_norm_dt6(10) ];

% At t = 8
dt_array = [ 0.01, 0.05, 0.1, 0.25, 0.5, 1.0 ];
L2_norm_final_T = [ L2_norm_dt1(751), L2_norm_dt4(151), ...
    L2_norm_dt2(76), L2_norm_dt5(31), L2_norm_dt3(16), L2_norm_dt6(8) ];

loglog( dt_array, L2_norm_final_T, 'o-', 'Color', red, ...
    'LineWidth', 2, ...
    'MarkerFaceColor', red, 'MarkerSize', 8 )
hold on
loglog( dt_array, 5e-2 * dt_array.^1, '--', 'Color', gray2 )
grid on
legend( 'error', '$ \mathcal{O}(h_{t}) $', 'FontSize', 14, 'Location', 'SE' )
xlabel('$ h_{t} $')
ylabel('$ \| w_{T} - w_{\mathrm{ref},T} \|_{L^{2}(\Omega)} $')

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                         Save figure 2                        |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = 'plots/ACE_stability_barrier';

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName_plot);
%--------------------------------------------------------------------------
