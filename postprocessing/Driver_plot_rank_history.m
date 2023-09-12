%==========================================================================
% Driver for plotting the

% Created:     2023.01.30
% Last change: 2023.02.13

%   Feb 13, 2023:
%       Added code for the plotting the full rank history during the first
%       few seconds.
%   Jan 30, 2023:
%       Created.
%==========================================================================

close all; clear; clc;

options_plot;

%--------------------------------------------------------------------------

% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K11_dt0.05.mat")

full_rank = size(W_hist_lra(1).U,1);
dt = t_hist_lra(2) - t_hist_lra(1);

t_hist_full_rank = 0:dt:t_hist_lra(1)-dt;
rank_W_hist_full_rank = full_rank * ones(1,length(t_hist_full_rank));

t_hist_lra_combo = [ t_hist_full_rank, t_hist_lra ];
rank_W_hist_lra_combo = [ rank_W_hist_full_rank, rank_W_hist_lra ];

%--------------------------------------------------------------------------
% Plot
semilogy( t_hist_lra_combo, rank_W_hist_lra_combo, "Color", red, "LineStyle", "-", "LineWidth", 2 )
grid on
xlabel('$t$')
ylabel('$ \mathrm{rank}(W_{\mathrm{ref}}) $')
% xlim( [ 0, 3 ] )
ylim( [ min(rank_W_hist_lra)-3  full_rank+10 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save figure                        |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = [ '../plots/ACE_lra_256x256_rank_vs_time_dt', num2str(dt) ];

set(gcf, 'Color', 'w');
export_fig(fileName_plot, '-pdf', '-opengl','-r600');

% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
% fprintf('Saved graph to file %s.eps.\n', fileName_plot);
%--------------------------------------------------------------------------