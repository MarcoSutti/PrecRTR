%==========================================================================
% Driver for plotting the initial conditions and the final numerical
% solutions for 1000 realizations of the stochastic Fisher-KPP equation.
% This is a reproduction of the results in the preprint:
% [1] Stable rank-adaptive dynamically orthogonal Runge-Kutta schemes,
%     Charous, A. and Lermusiaux, P., available on arXiv, 2022.

% Created:     2023.03.22
% Last change: 2023.05.09

%   Mar 22, 2023:
%       Created.
%==========================================================================

close all; clear; clc;

options_plot;

% Load the saved data. Initial conditions and solution at final time T.
load("FKPP_ICs_Nx1000_Nr1000.mat")
load("FKPP_CNLF_W_hist_Nx1000_Nr1000_T10_Nt1601.mat")

WT = W_CNLF_hist_stride(:,:,end);

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% LineWidth for plot:
linewidth = 1;

% Interval boundaries:
wL = 0;
wR = 40;
%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------

% Total number of realizations:
Nr = size(W0, 2);

% Spatial discretization:
Nx = size(W0, 1);   % number of spatial grid points

x = linspace( wL, wR, Nx)';

T = 10;

% Random colors in the violet-indigo palette:
% rnd_color_violet_indigo = [ rand(Nr,1), zeros(Nr,1), (rand/2+0.5) * ones(Nr,1) ];
% [0.667, 0.706, 0.118];

% Random colors in the yellow-orange palette:
rnd_color_green_brown = [ 0.667*ones(Nr,1), rand(Nr,1), 0.118*rand(Nr,1) ];

% rnd_color_yellow_orange = [ ones(Nr,1), rand(Nr,1), 0.259*ones(Nr,1) ];

%--------------------------------------------------------------------------
% Plot the initial condition:
figure(1)
plot( x, W0, "LineWidth", linewidth )
colororder( rnd_color_green_brown )
grid on
xlabel('$ x $')
ylabel('$ w $')
xlim( [ wL 25 ] )
ylim( [ -0.02 1.02 ] )
%--------------------------------------------------------------------------
% SAVE IMAGE TO PDF FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                          Save figure                         |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = [ '../plots/FKPP_CNLF_W0_Nx', num2str(Nx), '_Nr', ...
    num2str(Nr), '_T', num2str(T) ];

set(gcf, 'Color', 'w');
export_fig(fileName_plot, '-pdf' ); %, '-opengl','-r600');

% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.pdf.\n', fileName_plot);
%--------------------------------------------------------------------------
pause(2)
%--------------------------------------------------------------------------
% Plot numerical solution at final time T:
figure(2)
plot( x, WT, "LineWidth", linewidth )
colororder( rnd_color_green_brown )
% colororder( [ rnd_color_violet_indigo; rnd_color_green_brown ] )
grid on
xlabel('$ x $')
ylabel('$ w $')
xlim( [ wL 25 ] )
ylim( [ -0.02 1.02 ] )
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
fileName_plot = [ '../plots/FKPP_CNLF_WT_Nx', num2str(Nx), '_Nr', ...
    num2str(Nr), '_T', num2str(T) ];

set(gcf, 'Color', 'w');
export_fig(fileName_plot, '-pdf' ); %, '-opengl','-r600');

% % Save plot to eps file
% saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.pdf.\n', fileName_plot);
%--------------------------------------------------------------------------
