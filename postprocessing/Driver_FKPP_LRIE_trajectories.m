%==========================================================================
% Driver for plotting the low-rank time evolution of the trajectories
% of the Fisher-KPP equation generated by 'Driver_PrecRTR_Fisher_KPP'.

% Created:     2023.04.06
% Last change: 2023.04.19

%   Apr 6, 2023:
%       Created.
%==========================================================================

close all; clear; clc;

options_plot;

% Load the saved data:
load("../mdatafiles/FKPP_Prec1_W_hist_lra_100x100_K15_dt0.1.mat")

%--------------------------------------------------------------------------
% Spatial discretization:
% Interval boundaries:
wL = 0;
wR = 40;
%--------------------------------------------------------------------------

% Spatial discretization:
Nx = size(W_hist_lra(1).U, 1);

x = linspace( wL, wR, Nx)';

% Colors in the yellow-orange palette
rnd_color_yellow_orange = [ ones(Nx,1), rand(Nx,1), 0.259*ones(Nx,1) ];

dt = t_hist_lra(2) - t_hist_lra(1);
%--------------------------------------------------------------------------

for iter=1:size(W_hist_lra,2)

    current_t = (iter-1)*dt;

    % Form the dense matrix:
    W_mat = W_hist_lra(iter).U * W_hist_lra(iter).S * W_hist_lra(iter).V';

    % Plot numerical solution at current time:
    plot( x, W_mat, "LineWidth", 1 )
    grid on
    % Random colors in the yellow-orange palette:
    colororder( rnd_color_yellow_orange )
    xlabel('$ x $')
    ylabel('$ w $')
    xlim( [ wL 25 ] )
    ylim( [ -0.02 1.02 ] )
%     title(['Low-rank+IE, rank = ', num2str(size(W_hist_lra(iter).S,1)), ...
%         ', iter ',num2str(iter),', time = ', num2str(current_t)])
    drawnow;

end

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                          Save figure                         |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = '../plots/FKPP_LRIE_all_trajectories_at_final_time';

set(gcf, 'Color', 'w');
export_fig(fileName_plot, '-pdf' ); %, '-opengl','-r600');

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.pdf.\n', fileName_plot);