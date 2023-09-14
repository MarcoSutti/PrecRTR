%==========================================================================
% Driver for plotting the time evolution of 
%    \| W - W_{\mathrm{ref}} \|_{L^{2}(\Omega)}
% for the Fisher-KPP equation, for varying time-step sizes.
% h = 0.025, 0.0125, 0.00625.
% Nt = 401, 801, 1601.

% Created:     2023.05.08
% Last change: 2023.09.14

%   May 8, 2023:
%       Created.
%==========================================================================

close all; clear; clc;

options_plot;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Spatial discretization:
% Interval boundaries:
wL = 0;
wR = 40;

% Spatial discretization:
Nx = 1000;   % number of spatial grid points

x = linspace( wL, wR, Nx)';

% Mesh size
hx = (wR-wL)/(Nx-1);
%--------------------------------------------------------------------------
% Time discretization
T = 10;   % final time
%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------

% Load the saved data of the IMEX-CNLF solution:
% load('FKPP_CNLF_W_hist_Nx1000_Nr1000_T10_Nt3201.mat')

%--------------------------------------------------------------------------
% Load the saved data of the LR-CNLF solution:
% Nt = 401
load("FKPP_RTR_W_hist_Nx1000_Nr1000_T10_Nt401_final_rank12.mat")
load("FKPP_CNLF_W_hist_Nx1000_Nr1000_T10_Nt401.mat")

Nt401 = size(rank_W_RTR_hist, 2);
stride_Nt401 = floor(Nt401/100);

Nt_reduced = ceil(Nt401/stride_Nt401);

% Initializations.
L2_norm_Nt401 = zeros( 1, Nt_reduced );

k = 0;
for i = 1:stride_Nt401:Nt401
    k = k + 1;

    % Dummy variable:
    Wn_RTR_hist_i_mat = Wn_RTR_hist(i).U * Wn_RTR_hist(i).S * Wn_RTR_hist(i).V';

    % Calculate L2-norms:
    L2_norm_Nt401(k) = norm( W_CNLF_hist_stride(:,:,k) - Wn_RTR_hist_i_mat, 'fro' ) * sqrt(hx);

end

%--------------------------------------------------------------------------
% Load the saved data of the LR-CNLF solution:
% Nt = 801
load("FKPP_RTR_W_hist_Nx1000_Nr1000_T10_Nt801_final_rank12.mat")
load("FKPP_CNLF_W_hist_Nx1000_Nr1000_T10_Nt801.mat")

Nt801 = size(rank_W_RTR_hist, 2);
stride_Nt801 = floor(Nt801/100);

Nt_reduced = ceil(Nt801/stride_Nt801);

% Initializations.
L2_norm_Nt801 = zeros( 1, Nt_reduced );

k = 0;
for i = 1:stride_Nt801:Nt801
    k = k + 1;

    % Dummy variable:
    Wn_RTR_hist_i_mat = Wn_RTR_hist(i).U * Wn_RTR_hist(i).S * Wn_RTR_hist(i).V';

    % Calculate L2-norms:
    L2_norm_Nt801(k) = norm( W_CNLF_hist_stride(:,:,k) - Wn_RTR_hist_i_mat, 'fro' ) * sqrt(hx);

end

%--------------------------------------------------------------------------
% Load the saved data of the LR-CNLF solution:
% Nt = 1601
load("FKPP_RTR_W_hist_Nx1000_Nr1000_T10_Nt1601_final_rank12.mat")
load("FKPP_CNLF_W_hist_Nx1000_Nr1000_T10_Nt1601.mat")

Nt1601 = size(rank_W_RTR_hist, 2);
stride_Nt1601 = floor(Nt1601/100);

Nt_reduced = ceil(Nt1601/stride_Nt1601);

% Initializations.
L2_norm_Nt1601 = zeros( 1, Nt_reduced );

k = 0;
for i = 1:stride_Nt1601:Nt1601
    k = k + 1;

    % Dummy variable:
    Wn_RTR_hist_i_mat = Wn_RTR_hist(i).U * Wn_RTR_hist(i).S * Wn_RTR_hist(i).V';

    % Calculate L2-norms:
    L2_norm_Nt1601(k) = norm( W_CNLF_hist_stride(:,:,k) - Wn_RTR_hist_i_mat, 'fro' ) * sqrt(hx);

end

%--------------------------------------------------------------------------
%% Plot
figure(1)
close all

semilogy( t_hist_stride, L2_norm_Nt401, "Color", red, "LineStyle", "--", "LineWidth", 3 )
grid on
hold on
semilogy( t_hist_stride, L2_norm_Nt801, "Color", blue, "LineStyle", ":", "LineWidth", 3 )
semilogy( t_hist_stride, L2_norm_Nt1601, "Color", green, "LineStyle", "-", "LineWidth", 3 )
%--------------------------------------------------------------------------
% MS, 2023.03.06:
% Change the color of the grid:
% set( gca,'GridColor', super_light_gray )

% Added on 2023.03.06:
% Questo è il background "esterno" alla figura
% set( gca,'Color', gray3 )
%--------------------------------------------------------------------------
xlabel('$t$')
ylabel('$ \| w_{\mathrm{RTR}}  - w_{\mathrm{CNLF}} \|_{L^{2}(\Omega)} $')
legend( '$ h = 0.025 $', ...
    '$ h = 0.0125 $', ...
    '$ h = 0.00625 $', ...
    'FontSize', 14, ...
    'Location', 'SE' )
xlim( [ t_hist_stride(1), t_hist_stride(end) ] )
ylim( [ 8e-10 3e-7 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                          Save figure                         |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = 'plots/FKPP_L2_norm_W_RTR_vs_W_CNLF_varying_dt';

set(gcf, 'Color', 'w');
export_fig(fileName_plot, '-pdf', '-opengl','-r600');

% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.pdf.\n', fileName_plot);
