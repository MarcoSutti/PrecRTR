%==========================================================================
% Driver for plotting results related to PrecRTR.

% Generate Figure 1 of
%       An efficient preconditioner for the Riemannian trust-regions
%       method on the manifold of fixed-rank matrices,
%       M. Sutti, B. Vandereycken, Tech. report (submitted), 2022.
%       https://arxiv.org/abs/....

% Created:     29.08.2022
% Last change: 30.08.2022

%   Aug 29, 2022:
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
problem_type = 'LYAP';   % problem type
lev = 12;                % finest level
K = 5;
precon = [0, 1];            % Preconditiong: 0 for No, 1 for Yes.
%--------------------------------------------------------------------------
% Convergence of err-W for different ranks
figure(1);

% Line colors
clr(1, :) = red; 
clr(2, :) = green;

% Line style
style = {'^--', 'o-'};

% Stride for controlling the marker and ticks position frequency; i.e.,
% plot a marker and a tick every stride points.
stride = 5;

fprintf('+--------------------------------------------------------------+\n');
fprintf('|                            Plotting                          |\n');
fprintf('+--------------------------------------------------------------+\n');
for i = precon
    fileName = ['../mdatafiles/', problem_type, '_Prec', num2str(i), ...
        '_Grad_', num2str(lev), '_K', num2str(K), '.mat'];
    
    load( fileName );
        
    normalized_gradient = tr_gradnorm/tr_gradnorm(1);
%     fprintf('+--------------------------------------------------------------+\n');
%     fprintf( 'lev = %d, R-grad(end) = %.8e\n', levels(i), normalized_gradient(end) )
        
    % Plot convergence behavior of gradient norm
    handle_array(i+1) = semilogy( 1:ntot_iter+1, normalized_gradient, style{i+1}, ...
        'Color', clr(i+1, :), 'LineWidth', 1.5, ...
        'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', clr(i+1, :), ...
        'MarkerSize', 8, ...
        'MarkerIndices', 1:stride:ntot_iter+1 );
    
    grid on
    hold on
    
end

% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;

drawnow;
for i=1:2
    handle_array(i).MarkerHandle.LineWidth = myMarkerLineWidth;
end

% Legend
handleLegend = legend( handle_array, ...
    {'without Prec', ...
     'with Prec'}, ...
    'FontSize', 11, 'Location', 'southeast' );

drawnow;
for i=1:2
    lineEntry = findobj(handleLegend.EntryContainer, 'Object', handle_array(i) );
    entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    entryMarker.LineWidth = myMarkerLineWidth;
end

% Axes labels and title
xlabel('iteration $ i $ of RTR');
ylabel('R-grad')
% ylim( [ 1e-15, 10 ] )

% xticks( 0:25:length() );
% xticklabels( 0:25:length() );

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save plot                          |\n');
fprintf('+--------------------------------------------------------------+\n');

fileName_plot = [ '../plots/', problem_type, '_Comparison_PrecRTR_', ...
    num2str(lev), '_K', num2str(K) ];

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName_plot);