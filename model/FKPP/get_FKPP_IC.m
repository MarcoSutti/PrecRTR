function [ W0, Romega ] = get_FKPP_IC( Nx, Nr, x )

% function [ W0, Romega ] = get_FKPP_IC( Nx, Nr, x )
% Purpose: Computes the initial condition for the Fisher-KPP equation.

% Created:     2023.05.09
% Last change: 2023.05.09

%   May 9, 2023:
%       Created.

% Initializations:
W0 = zeros(Nx, Nr);
r_omega_array = zeros(Nr, 1);  % Legacy: r_omega_array = zeros(1, tot_num_realizations);

for idx_realization=1:Nr
    %     fprintf( "Realization: %.d \n", idx_realization );

    % The coefficients a, b, r follow a uniform law.
    a_omega = unifrnd(1/5, 2/5);
    b_omega = unifrnd(1/10, 11/10);
    r_omega = unifrnd(1/4, 1/2);   % reaction rate

    % 1.1. Initial condition:
    w0 = a_omega * exp(-b_omega * x.^2);

    % Store w0 in W_CNLF_hist and W_EuTR_hist for later saving.
    %     time_iter = 1;

    W0(:,idx_realization) = w0;
    %     W_EuTR_hist(:,idx_realization,time_iter) = w0;

    % Store r_omega in r_omega_array for later saving
    r_omega_array(idx_realization) = r_omega;
end

Romega = spdiags(r_omega_array, 0, Nx, Nr);

end