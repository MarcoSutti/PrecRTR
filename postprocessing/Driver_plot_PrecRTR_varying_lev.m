%==========================================================================
% Driver for plotting PrecRTR varying the mesh-size.

% Generate Figure 1 of
%       An efficient preconditioner for the Riemannian trust-regions
%       method on the manifold of fixed-rank matrices,
%       M. Sutti, B. Vandereycken, Tech. report (submitted), 2022.
%       https://arxiv.org/abs/....

% Created:     30.08.2022
% Last change: 30.08.2022

%   Aug 30, 2022:
%       Created.
%==========================================================================
close all; clear; clc;

% startup;

plot_startup;

% Get the present working directory
str = pwd;

if contains(str,'experiments') == 0
    error('Please make sure that you are running this script in the ''PrecRTR/experiments'' folder.')
end
%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
problem_type = 'NPDE';   % problem type
levels = 10:15;          % levels of discretization
K = 5;
precon = 1;              % Preconditiong: 0 for No, 1 for Yes.
%--------------------------------------------------------------------------
% Convergence of err-W for different ranks
figure(1);

% Line colors
clr(1, :) = red; 
clr(2, :) = green;
clr(3, :) = blue;
clr(4, :) = yellow;
clr(5, :) = orange;
clr(6, :) = darkgray;

% Line style
style = {'^--', 'o-'};

% Stride for controlling the marker and ticks position frequency; i.e.,
% plot a marker and a tick every stride points.
stride = 5;

fprintf('+--------------------------------------------------------------+\n');
fprintf('|                            Plotting                          |\n');
fprintf('+--------------------------------------------------------------+\n');
for i = 1:length(levels)
    fileName = ['../mdatafiles/', problem_type, '_Prec', num2str(precon), ...
        '_Grad_', num2str(levels(i)), '_K', num2str(K), '.mat'];
    
    load( fileName );
    
    normalized_gradient = tr_gradnorm/tr_gradnorm(1);
    n_outer = length(normalized_gradient);
    
    fprintf('+--------------------------------------------------------------+\n');
    fprintf( 'lev = %d, R-grad(end) = %.8e\n', levels(i), normalized_gradient(end) )
    
    % Plot convergence behavior of gradient norm
    handle_array(i) = semilogy( 1:n_outer, normalized_gradient, 'o-', ...
        'Color', clr(i, :), 'LineWidth', 2, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', clr(i, :), ...
        'MarkerSize', 8, ...
        'MarkerIndices', 1:stride:n_outer );
    
    grid on
    hold on
    
end

% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;

drawnow;
for i=1:length(levels)
    handle_array(i).MarkerHandle.LineWidth = myMarkerLineWidth;
end

% Legend
handleLegend = legend({num2str(levels(1)), num2str(levels(2)), ...
    num2str(levels(3)), num2str(levels(4)), num2str(levels(5)), num2str(levels(6))}, ...
    'Location','SW');
handleLegend.Title.String = '$\ell_{\mathrm{f}}$';

drawnow;
for i=1:length(levels)
    lineEntry = findobj(handleLegend.EntryContainer, 'Object', handle_array(i) );
    entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    entryMarker.LineWidth = myMarkerLineWidth;
end


% Axes labels and title
xlabel('iteration $ i $ of RTR');
ylabel('R-grad')
ylim( [ 1e-10, 1e3 ] )

% xticks( 0:25:length() );
% xticklabels( 0:25:length() );

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save plot                          |\n');
fprintf('+--------------------------------------------------------------+\n');

fileName_plot = [ '../plots/', problem_type, '_PrecRTR_level_', ...
    num2str(levels(1)), '_to_', num2str(levels(end)), '_K', num2str(K) ];

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName_plot);
