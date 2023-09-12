function [ t, stat_hzls ] = getLinesearch(problem, x, d, storedb, key)

% Marco: I added the output structure stat_hzls to keep track of the number
% of cost function evaluations.

% Returns a hint for line-search algorithms.
%
% function t = getLinesearch(problem, x, d)
% function t = getLinesearch(problem, x, d, storedb)
% function t = getLinesearch(problem, x, d, storedb, key)
%
% For a line-search problem at x along the tangent direction d, computes
% and returns t such that retracting t*d at x yields a good point around
% where to look for a line-search solution. That is: t is a hint as to
% "how far to look" along the line.
% 
% storedb is a StoreDB object, key is the StoreDB key to point x.
%
% See also: canGetLinesearch

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, July 17, 2014.
% Contributors: 
% Change log: 
%
%   April 3, 2015 (NB):
%       Works with the new StoreDB class system.

    % Allow omission of the key, and even of storedb.
    if ~exist('key', 'var')
        if ~exist('storedb', 'var')
            storedb = StoreDB();
        end
        key = storedb.getNewKey();
    end


    if isfield(problem, 'linesearch')
    %% Compute the line-search hint function using linesearch.
    
        % Check whether this function wants to deal with storedb or not.
        switch nargin(problem.linesearch)
            case 2
                t = problem.linesearch(x, d);
            case 3
                % Obtain, pass along, and save the store for x.
                store = storedb.getWithShared(key);
                [ t, stat_hzls ] = problem.linesearch(problem, x, d);
                % Marco, 15.11.2019: I added the output structure stat_hzls to get
                % the number of function evaluations of Hager-Zhang line search out
                % of the function.
                % Original version: [ t ] = problem.linesearch(problem, x, d);
                storedb.setWithShared(store, key);
            case 4
                % Pass along the whole storedb (by reference), with key.
                t = problem.linesearch(x, d, storedb, key);
            otherwise
                up = MException('manopt:getLinesearch:badfun', ...
                    'linesearch should accept 2, 3 or 4 inputs.');
                throw(up);
        end

    else
    %% Abandon computing the line-search function.

        up = MException('manopt:getLinesearch:fail', ...
            ['The problem description is not explicit enough to ' ...
             'compute a line-search hint.']);
        throw(up);
        
    end
    
end
