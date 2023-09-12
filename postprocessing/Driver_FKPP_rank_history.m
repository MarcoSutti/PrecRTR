%==========================================================================
% Driver for plotting the rank history of the LRIE solution to the
% Fisher-KPP equation.

% Created:     2023.04.10
% Last change: 2023.04.11

%   Apr 10, 2023:
%       Created.
%==========================================================================

close all; clear; clc;

options_plot;

%--------------------------------------------------------------------------

% Load the saved data:
% load("FKPP_Prec0_W_hist_lra_1000x1000_K9_dt0.125.mat")
% load("FKPP_Prec0_W_hist_lra_1000x1000_K22_dt0.125.mat")
load("../mdatafiles/FKPP_Prec1_W_hist_lra_100x100_K15_dt0.1.mat")

dt = t_hist_lra(2) - t_hist_lra(1);

%--------------------------------------------------------------------------
% Plot
plot( t_hist_lra, rank_W_hist_lra, "Color", red, "LineStyle", "-", "LineWidth", 2 )
grid on
xlabel('$t$')
ylabel('$ \mathrm{rank}(W) $')
xlim( [ t_hist_lra(1), t_hist_lra(end) ] )
ylim( [ 6  max(rank_W_hist_lra)+1 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save figure                        |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = [ '../plots/FKPP_lra_rank_vs_time_dt', num2str(dt) ];

set(gcf, 'Color', 'w');
export_fig(fileName_plot, '-pdf', '-opengl','-r600');

% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.pdf.\n', fileName_plot);
%--------------------------------------------------------------------------