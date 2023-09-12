%==========================================================================
% Driver for plotting PrecRTR varying low rank.
% Driver for plotting on the convergence behavior of the quantity err-W
% % Approximate solution Wh: relative error between a reference solution 
% Wh_ref and the approximate solution Wh given by RMGLS.

% Generate Figure --- of
%       An efficient preconditioner for the Riemannian trust-regions
%       method on the manifold of fixed-rank matrices,
%       M. Sutti, B. Vandereycken, Tech. report (submitted), 2022.
%       https://arxiv.org/abs/....

% Created:     2022.09.28
% Last change: 2022.09.28

%   Sept 28, 2022:
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
lev = 10;                % level of discretization
ranks = 2:2:10;
precon = 1;              % Preconditiong: 0 for No, 1 for Yes.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% fprintf('+--------------------------------------------------------------+\n');
% fprintf('|          Load reference solution (Euclidean sol.)...         |\n');
% fprintf('+--------------------------------------------------------------+\n');
% 
% [ Wh_ref, norm_Wh_ref ] = get_Wh_ref( pars );
% 
% Wh_0_mat = Wh.U * Wh.S * Wh.V';
% 
% % err-W at the 0th iteration:
% err_W_0 = norm( Wh_0_mat - Wh_ref, 'fro' )/norm_Wh_ref;
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
for i = 1:length(ranks)
    fileName = ['../mdatafiles/', problem_type, '_Prec', num2str(precon), ...
        '_Grad_', num2str(lev), '_K', num2str(ranks(i)), '.mat'];
    
    load( fileName );
    
    normalized_gradient = tr_gradnorm/tr_gradnorm(1);
    n_outer = length(normalized_gradient);
    
    fprintf('+--------------------------------------------------------------+\n');
    fprintf( 'rank = %d, R-grad(end) = %.8e\n', ranks(i), normalized_gradient(end) )
    
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
for i=1:length(ranks)
    handle_array(i).MarkerHandle.LineWidth = myMarkerLineWidth;
end

% Legend
handleLegend = legend({num2str(ranks(1)), num2str(ranks(2)), ...
    num2str(ranks(3)), num2str(ranks(4)), num2str(ranks(5)) }, ...
    'Location','SW');
handleLegend.Title.String = '$ k $';

drawnow;
for i=1:length(ranks)
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

fileName_plot = [ '../plots/', problem_type, '_PrecRTR_rank_', ...
    num2str(ranks(1)), '_to_', num2str(ranks(end)), '_lev', num2str(lev) ];

% Save plot to eps file
saveas( gcf, fileName_plot, 'epsc' )
fprintf('Saved graph to file %s.eps.\n', fileName_plot);
