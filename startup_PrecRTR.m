%--------------------------------------------------------------------------
% Startup PrecRTR
% Created:     2022.08.29
% Last change: 2022.09.29
%--------------------------------------------------------------------------
close all; clear; clc;

folder_RMGLS = '../../Marco_Multilevel/matlab/RMGLS-release01/1_RMGLS_EML';
folder_PrecRTR = '../PrecRTR_matlab';

addpath( genpath(folder_RMGLS) );
addpath( genpath(folder_PrecRTR) );

% MS, 30.05.2020: The following piece of code is taken over from the manopt
%                 script 'importmanopt'
% Ask user if the path should be saved or not
fprintf('1_RMGLS_EML and PrecRTR were added to Matlab''s path.\n');
response = input('Save path for future Matlab sessions? [Y/N] ', 's');
if strcmpi(response, 'Y')
    failed = savepath();
    if ~failed
        fprintf('Path saved: no need to call startup next time.\n');
    else
        fprintf(['Something went wrong.. Perhaps missing permission ' ...
                 'to write on pathdef.m?\nPath not saved: ' ...
                 'please re-call startup_RMGLS next time.\n']);
    end
else
    fprintf('Path not saved: please re-call startup_PrecRTR next time.\n');
end


% Install manopt:
cd manopt;
importmanopt;
cd ..;
