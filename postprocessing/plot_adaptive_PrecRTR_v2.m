function [ ] = plot_adaptive_PrecRTR_v2( res_W_adaptive, res_W_no_adapt, err_W_adaptive, ...
    err_W_no_adapt, idmg_change_rank, pars )

% function [ ] = plot_adaptive_PrecRTR_v2( res_W_adaptive, res_W_no_adapt, err_W_adaptive, ...
%     err_W_no_adapt, idmg_change_rank, pars )
% Purpose: Plot convergence behavior for the adaptive PrecRTR, v2.

% Created:     2022.10.18
% Last change: 2022.10.18

%   Oct 18, 2022:
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
handle_array(5) = semilogy( iters_adaptive( idmg_change_rank == 1 ), res_W_adaptive( idmg_change_rank == 1 ), ...
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
handle_array(6) = semilogy( iters_adaptive( idmg_change_rank == 1 ), err_W_adaptive( idmg_change_rank == 1 ), ...
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
    'k = 10, $ r(W) $', ...
    'k = 10, err-$ W $', ...
    'Location','SW', 'FontSize', 12);
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
ylim( [ 1e-13, 1e2 ] )

%--------------------------------------------------------------------------
% SAVE IMAGE TO EPS FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save plot                          |\n');
fprintf('+--------------------------------------------------------------+\n');

fileName_plot = [ 'plots/', pars.problem_type, '_Adaptive_v2_rank_5_to_10_lev', num2str(pars.lev_finest) ];

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName_plot);

end