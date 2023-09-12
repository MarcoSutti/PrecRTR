function [ ] = plot_adaptive_PrecRTR( res_W_adaptive, res_W_no_adapt, ...
    err_W_adaptive, err_W_no_adapt, idmg_red, ranks_array, pars )

% function [ ] = plot_adaptive_PrecRTR( res_W_adaptive, res_W_no_adapt, ...
%    err_W_adaptive, err_W_no_adapt, idmg_red, ranks_array, pars )
% Purpose: Plot convergence behavior for the adaptive PrecRTR.
% Created:     2022.09.29
% Last change: 2022.10.18

%   Oct 18, 2022:
%       Minor changes to the format of the legend.
%   Sept 29, 2022:
%       Created.

plot_startup;

LineWidth = 2.0;
MarkerSizeChangeInRank = 10;
MarkerColorChangeInRank = yellow;

% % LineWidth of the MarkerEdge:
% myMarkerLineWidth = 0.5;

fprintf('+--------------------------------------------------------------+\n');
fprintf('|                            Plotting                          |\n');
fprintf('+--------------------------------------------------------------+\n');

% 1) Adaptive PrecRTR
% Plot convergence behavior of residual norm
iters_adaptive = 1:length(res_W_adaptive);
iters_no_adapt = 1:length(res_W_no_adapt);

handle_array(1) = semilogy( iters_adaptive, res_W_adaptive, '-', ...
    'Color', Cherry, ...
    'LineWidth', LineWidth );

grid on
hold on

% Convergence behavior of err-W
handle_array(2) = semilogy( iters_adaptive, err_W_adaptive, '--', ...
    'Color', Cherry, ...
    'LineWidth', LineWidth );

% Mark the points where the rank changes, for res_W_adaptive.
handle_array(5) = semilogy( iters_adaptive( idmg_red == 1 ), res_W_adaptive( idmg_red == 1 ), ...
    's', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', MarkerColorChangeInRank, ...
    'MarkerSize', MarkerSizeChangeInRank );

% 2) Plain (i.e., non-adaptive) PrecRTR
handle_array(4) = semilogy( iters_no_adapt, res_W_no_adapt, '-', ...
    'Color', Teal, ...
    'LineWidth', LineWidth );

handle_array(5) = semilogy( iters_no_adapt, err_W_no_adapt, '--', ...
    'Color', Teal, ...
    'LineWidth', LineWidth );

% Mark the points where the rank changes, for err_W_adaptive.
handle_array(6) = semilogy( iters_adaptive( idmg_red == 1 ), err_W_adaptive( idmg_red == 1 ), ...
    's', ...
    'MarkerEdgeColor', 'k', ...
    'MarkerFaceColor', MarkerColorChangeInRank, ...
    'MarkerSize', MarkerSizeChangeInRank );

%
% drawnow;
% for i=1:4
%     handle_array(i).MarkerHandle.LineWidth = myMarkerLineWidth;
% end

% Legend
handleLegend = legend('Adaptive, $ r(W) $', ...
    'Adaptive, err-$ W $', ...
    'Change in rank', ...
    ['k = ', num2str(pars.K) ', $ r(W) $'], ...
    ['k = ', num2str(pars.K) ', err-$ W $'], ...
    'Location', 'SE', 'FontSize', 12 );
% handleLegend.Title.String = 'Method';

% drawnow;
% for i=1:4
%     lineEntry = findobj(handleLegend.EntryContainer, 'Object', handle_array(i) );
%     entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
%     entryMarker.LineWidth = myMarkerLineWidth;
% end

% Axes labels and title
xlabel('iteration $ i $ of PrecRTR');
ylabel('$ r(W) $ and err-$ W $')
xlim( [ 0, iters_adaptive(end) ] )
ylim( [ 1e-15, 1e1 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save plot                          |\n');
fprintf('+--------------------------------------------------------------+\n');

fileName_plot = [ 'plots/', pars.problem_type, '_Adaptive_rank_', ...
    num2str(ranks_array(1)), '_to_', num2str(ranks_array(end)), '_lev', num2str(pars.lev_finest) ];

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName_plot);

end