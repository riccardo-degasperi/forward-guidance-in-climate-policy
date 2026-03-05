%-------------------------------------------------------------------------%
% Define sign restrictions
%
% The SR matrix is a 3-D matrix with dimension (nvar,3,nshocks). Assume you
% have 3 variables and you want to identify just one shock. Then:
% 
%              from        to          sign
%  SR(:,:,1) = [ 1           4           1          % VAR1
%                1           4          -1          % VAR2
%                0           0           0];        % VAR3
% 
%   - Column 1 defines the 1st period from which the restriction is imposed 
%   - Column 2 defines the last period to which the restriction is imposed 
%   - Column 3 defines the sign of the restriction: +ve (1) or -ve (-1)
%   - To leave unrestricted set all three columns to zero
% 
% In the above example we set the following restrictions. VAR1 must respond
% with positive sign from period 1 to period 4; VAR2 must respond with 
% negative sign from period 1 to period 4; VAR3 is left unrestricted
%
% An additional shock could be defined as follows:
%              from        to          sign
%  SR(:,:,2) = [ 1           4           1          % VAR1
%                1           4           1          % VAR2
%                1           4          -1];        % VAR3
%
% If you identify some shocks with proxies, it makes sense to impose
% additional sign restrictions on those shocks only if the parameters of
% the reduced form change at each iteration, in a full Gibbs sampler.
%
% Adapted from Ambrogio Cesa-Bianchi, March 2015
%
% r.degasperi@warwick.ac.uk                                     22/09/2020
%-------------------------------------------------------------------------%


% Cost-push shock
SR(:,:,1) = [ 0     0     0    % EUA fq
              0     0     0    % Stoxx Europe 600
              0     0     0    % Critical Materials slope
              0     0     0    %
              0     0     0    %
              0     0     0    %
              0     0     0    %
              0     0     0    %
              1     1     1    %
              1     1    -1    %
              1     1    -1];  %

% Residual endogenous component
SR(:,:,2) = [ 0     0     0    % EUA fq
              0     0     0    % Stoxx Europe 600
              0     0     0    % Critical Materials slope
              0     0     0    %
              0     0     0    %
              0     0     0    %
              0     0     0    %
              0     0     0    %
              1     1     1    %
              1     1     1    %
              1     1     1];  %

% Path shock
SR(:,:,3) = [ 0     0     0    % EUA fq
              0     0     0    % Stoxx Europe 600
              0     0     0    % Critical Materials slope
              0     0     0    %
              0     0     0    %
              0     0     0    %
              0     0     0    %
              0     0     0    %
              1     1     1    %
              1     1    -1    %
              1     1     1];  %


%-------------------------------------------------------------------------%
% Pack output
modelOpt.SR = SR;
clear SR





