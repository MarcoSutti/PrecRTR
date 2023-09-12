function [ ] = plot_PrecRTR( tr_gradnorm, pars )

% function [ ] = plot_PrecRTR( tr_gradnorm, pars )

% Created:     2022.09.13
% Last change: 2022.09.13

%   Sept 13, 2022:
%       Created.

fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Plotting                           |\n');
fprintf('+--------------------------------------------------------------+\n');

plot_startup;

ntot_iter = length(tr_gradnorm);

% Stride for controlling the marker and ticks position frequency; i.e.,
% plot a marker and a tick every stride points.
stride = 5;

if pars.precon==1
    clr = red;
elseif pars.precon==0
    clr = green;
end

% Plot convergence behavior of gradient norm
h_semilogy = semilogy( 1:ntot_iter, tr_gradnorm, '^-', ...
    'Color', clr, 'LineWidth', 1.5, ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', clr, ...
    'MarkerSize', 8, ...
    'MarkerIndices', 1:stride:ntot_iter );

grid on


% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;

drawnow;
h_semilogy.MarkerHandle.LineWidth = myMarkerLineWidth;

% Legend
if pars.precon==0
    handleLegend = legend( h_semilogy, ...
        {'RTR'}, ...
        'FontSize', 11, 'Location', 'southeast' );
elseif pars.precon==1
    handleLegend = legend( h_semilogy, ...
        {'PrecRTR'}, ...
        'FontSize', 11, 'Location', 'southeast' );
end


drawnow;
lineEntry = findobj(handleLegend.EntryContainer, 'Object', h_semilogy );
entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
entryMarker.LineWidth = myMarkerLineWidth;

% Axes labels and title
xlabel('iteration $ i $ of RTR');
ylabel('R-grad')
ylim([eps, 1])


%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save plot                          |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = [ 'plots/', pars.problem_type, '_Prec', num2str(pars.precon), ...
    '_Grad_', num2str(pars.lev_finest), '_K', num2str(pars.K) ];

if pars.precon==0
    title('Without Hessian preconditioner')
elseif pars.precon==1
    title('With Hessian preconditioner')
end

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName_plot);

end
