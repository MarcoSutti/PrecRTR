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
load("ACE_Prec1_W_hist_lra_256x256_K11_dt0.5.mat")


%--------------------------------------------------------------------------
% Plot
semilogy( t_hist_lra, rank_W_hist_lra, "Color", 'r', "LineStyle", "-", "LineWidth", 2 )
hold on
grid on
xlabel('$t$')
ylabel('$ \mathrm{rank}(W_{\mathrm{ref}}) $')
% xlim( [ 0, 3 ] )
ylim( [ min(rank_W_hist_lra)-3 max(rank_W_hist_lra)+5 ] )