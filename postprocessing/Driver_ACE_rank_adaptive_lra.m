%==========================================================================
% Driver for plotting the numerical rank history of W for the Allen-Cahn
% equation from the rank-adaptive version.

% Created:     2023.01.30
% Last change: 2023.01.30

%   Jan 30, 2023:
%       Created.
%==========================================================================

close all; clear; clc;

options_plot;

% Load the saved data:
load("ACE_Prec1_W_hist_lra_256x256_K11_dt0.05.mat")


% Create bar for full-rank first 0.5 seconds:
t_hist_lra_ext = [ 0, 0.5, t_hist_lra ];
rank_W_hist_lra_ext = [ 256, 256, rank_W_hist_lra ];

%--------------------------------------------------------------------------
% Plot
semilogy( t_hist_lra_ext, rank_W_hist_lra_ext, "Color", red, ...
    "LineStyle", "-", "LineWidth", 2 )
hold on
grid on
xlabel('$t$')
ylabel('$ \mathrm{rank}(W_{\mathrm{ref}}) $')
xlim( [ 0, t_hist_lra(end) ] )
ylim( [ min(rank_W_hist_lra_ext)-3 max(rank_W_hist_lra_ext)+10 ] )


%--------------------------------------------------------------------------
% SAVE IMAGE TO PDF FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                          Save figure                         |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = [ '../plots/ACE_lra_256x256_rank_vs_time_dt0.05' ];

% ACE_lra_256x256_rank_vs_time_dt0.05.pdf

% Save plot to pdf file
export_fig(fileName_plot, '-pdf' );
fprintf('Saved graph to file %s.pdf.\n', fileName_plot);
%--------------------------------------------------------------------------