%==========================================================================
% Driver for the Allen-Cahn equation reference solution computed with the
% fourth-order explicit Runge-Kutta method (ERK4).

% Created:     2022.11.29
% Last change: 2023.02.10

%   Jan 10, 2023:
%       Added solution with the CFD discretization of the Laplacian.
%   Dec 31, 2022:
%       Added saving of the numerical solution every (ns * iter) seconds.
%   Dec 30, 2022:
%       Fixed mistake with periodic BCs.
%   Dec 29, 2022:
%       Created.
%==========================================================================

close all; clear; clc;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Spatial discretization:
Nx = 256; % Be aware of the stability condition for explicit methods...
Ny = Nx;
Lx = 2*pi;
Ly = 2*pi;
hx = Lx/Nx;
hy = Ly/Ny;
% Square domain as in the paper of Rodgers and Venturi.
x = linspace( 0, Lx, Nx+1);
y = linspace( 0, Ly, Ny+1);
x = x(1:end-1);
y = y(1:end-1);
[ xx , yy ] = ndgrid(x,y);

epsilon = 0.1; % value used in the paper of Rodgers and Venturi.

% Initial condition from the paper of Rodgers and Venturi, 2022:
W0 = get_ACE_IC_RV2022( xx, yy );

% Time discretization:
dt = 1e-4;
T = 0.5;
Nt = round( T/dt );

% Save the data every "ns" iterations:
ns = 1000;

% 2) The discretized 1D Laplacian in sparse format:
A = -(1/hx^2) * get_Ah( Nx );

plot = true;

%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------

fun = @(W) get_RHS( W, A, epsilon );

W = W0;

% MS, 2022.12.31: Save history of W every (ns * dt) seconds.
W_hist(:,:,1) = W0;
it_hist = 1;
t_hist = zeros( Nt/ns, 1 );

for iter=1:Nt

    k1 = fun( W );
    k2 = fun( W + (dt/2)*k1 );
    k3 = fun( W + (dt/2)*k2 );
    k4 = fun( W + dt*k3 );
    W = W + (dt/6)*( k1 + 2*k2 + 2*k3 + k4 );


    % Save history of the approximate solution every ns steps:
    if mod(iter,ns)==0
        it_hist = it_hist + 1;
        t_hist(it_hist) = iter*dt;
        W_hist(:,:,it_hist) = W;

        % For plotting purposes:
        if plot
            [ ~, handle_cf ] = contourf( x, y, W, 100 );
            set( handle_cf, 'edgecolor', 'none' );
            axis equal
            colormap("parula");
            colorbar;
            title(['Explicit RK4, iter ',num2str(iter),', time = ', num2str(dt*iter)])
            pause(0.01)
        end
    end
end

%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
%--------------------------------------------------------------------------
% SAVE DATA TO MAT-FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save data                          |\n');
fprintf('+--------------------------------------------------------------+\n');
% 2022.12.31: Save the history of t and W.
fileName_mfile = [ 'reference_solutions/ACE_ref_', num2str(Nx), 'x', num2str(Nx), '_T', ...
    num2str(T), '_dt', num2str(dt), '.mat'];

save( fileName_mfile, 't_hist', 'W_hist')
fprintf('Saved data to file %s.\n', fileName_mfile);
%--------------------------------------------------------------------------

RHS_at_T = get_RHS( W, A, epsilon );

fprintf('+--------------------------------------------------------------+\n');
fprintf( 'Residual at T = %.d : %.4e.\n', T, norm( RHS_at_T, 'fro' ) * hx );
