%==========================================================================
% Driver for plotting the time history of the rank of the reference
% solution the stochastic Fisher-KPP equation.
% The data are generated by 'Driver_compare_CNLF_vs_Eucl_TR_multiple_realizations.m'.

% Created:     2023.03.24
% Last change: 2023.05.09

%   Mar 24, 2023:
%       Created.
%==========================================================================

close all; clear; clc;

options_plot;

% Load the saved data generated by ''
load("FKPP_CNLF_W_hist_Nx1000_Nr1000_T10_Nt1601.mat")
load("FKPP_RTR_W_hist_Nx1000_Nr1000_T10_Nt1601_final_rank12.mat")

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Final time:
T = 10;

% Total number of time discretization points:
Nt = 1601;

% Total number of realizations:
Nr = 1000;

% Spatial discretization:
Nx = 1000;   % number of spatial grid points

t_hist = 0:T/(size(rank_W_RTR_hist,2)-1):T;
t_hist_stride = 0:T/(size(rank_W_CNLF_hist_stride,2)-1):T;
%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------

%%
close all; clc;

%--------------------------------------------------------------------------
% Plot rank evolution
figure(1)
plot( t_hist, rank_W_RTR_hist, "Color", Cerulean, "LineWidth", 3 )
hold on
plot( t_hist_stride, rank_W_CNLF_hist_stride, "--", "Color", green, "LineWidth", 3 )
grid on
xlabel('$ t $')
ylabel('numerical rank of $ W $')
legend( 'LR-CNLF', 'CNLF', 'FontSize', 14, 'Location', 'SE' )
xlim( [ t_hist(1)   t_hist(end) ] )
% ylim( [ rank_r-6   max(rank_hist_before_truncation)+2 ] )
ylim( [ min(rank_W_RTR_hist)-5   max(rank_W_CNLF_hist_stride)+2 ] )
% %--------------------------------------------------------------------------
% % Get rid of the large invisible border of the figure:
% set(gca, 'units', 'normalized'); % Just making sure it's normalized
% Tight = get(gca, 'TightInset');  % Gives you the bording spacing between plot box and any axis labels
% % [Left Bottom Right Top] spacing
% NewPos = [Tight(1) Tight(2) 1-Tight(1)-Tight(3) 1-Tight(2)-3*Tight(4)]; % New plot position [X Y W H]
% set(gca, 'Position', NewPos);
% %--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% SAVE IMAGE TO PDF FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                          Save figure                         |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = [ '../plots/FKPP_rank_history_Nx', num2str(Nx), '_Nr', ...
    num2str(Nr), '_T', num2str(T) ];

set(gcf, 'Color', 'w');
export_fig(fileName_plot, '-pdf', '-opengl','-r600');

% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.pdf.\n', fileName_plot);
%--------------------------------------------------------------------------
