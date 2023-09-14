%==========================================================================
% Driver for checking agreement between solutions computed on 2023.05.09
% and thos computed on September 2023.
% Last change: 2023.09.14

%   Sep 14, 2023:
%       Created.
%==========================================================================

close all; clear; clc;

% % Load ICs saved on 2023.05.09:
% load('FKPP_ICs_Nx1000_Nr1000_20230509.mat')
% 
% W0_20230509 = W0;
% Romega_20230509 = Romega;
% 
% % Load ICs saved on 2023.09.13:
% load('FKPP_ICs_Nx1000_Nr1000.mat')
% 
% % Check:
% norm( W0_20230509 - W0, 'fro')
% norm( Romega_20230509 - Romega, 'fro')



% load('FKPP_RTR_W_hist_Nx1000_Nr1000_T10_Nt1601_final_rank12_20230509.mat')
% rank_W_RTR_hist_20230509 = rank_W_RTR_hist;
% Wn_RTR_hist_20230509 = Wn_RTR_hist;
% 
% load('FKPP_RTR_W_hist_Nx1000_Nr1000_T10_Nt1601_final_rank12.mat')
% % Check:
% norm( rank_W_RTR_hist_20230509 - rank_W_RTR_hist, 'fro')
% 
% for i=1:length(rank_W_RTR_hist_20230509)
%     norm_U(i) = norm( Wn_RTR_hist_20230509(i).U - Wn_RTR_hist(i).U, 'fro');
%     norm_S(i) = norm( Wn_RTR_hist_20230509(i).S - Wn_RTR_hist(i).S, 'fro');
%     norm_V(i) = norm( Wn_RTR_hist_20230509(i).V - Wn_RTR_hist(i).V, 'fro');
% end
% 
% sum(norm_U)
% sum(norm_S)
% sum(norm_V)


% load('FKPP_RTR_W_hist_Nx1000_Nr1000_T10_Nt801_final_rank12_20230509.mat')
% load('FKPP_RTR_W_hist_Nx1000_Nr1000_T10_Nt1601_final_rank12_20230509.mat')


load('FKPP_CNLF_W_hist_Nx1000_Nr1000_T10_Nt3201_20230509.mat')

t_hist_stride_20230509 = t_hist_stride;
rank_W_CNLF_hist_stride_20230509 = rank_W_CNLF_hist_stride;
W_CNLF_hist_stride_20230509 = W_CNLF_hist_stride;

clear t_hist_stride
clear rank_W_CNLF_hist_stride
clear W_CNLF_hist_stride

load('FKPP_CNLF_W_hist_Nx1000_Nr1000_T10_Nt3201.mat')

% Check:
norm( rank_W_CNLF_hist_stride_20230509 - rank_W_CNLF_hist_stride, 'fro')
norm( t_hist_stride_20230509 - t_hist_stride, 'fro')

% In rank_W_CNLF_hist_stride_20230509 manca il rango della condizione iniziale.
% Lo aggiungo manualmente, e di conseguenza tutto sara' coerente nelle due
% versioni:
rank_W_CNLF_hist_stride_20230509(1) = 22;

for i=1:length(t_hist_stride_20230509)
    norm_diff_W(i) = norm( W_CNLF_hist_stride_20230509(i) - W_CNLF_hist_stride(i), 'fro');
end


sum(norm_diff_W)
