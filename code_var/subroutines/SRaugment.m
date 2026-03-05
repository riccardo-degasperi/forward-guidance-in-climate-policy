function SRnew = SRaugment(SR,dims,dataOpt)
% SRaugment checks that the sign restriction script has been compiled
% properly and prepares the restrictions for later use.
%
% last modified: 26/03/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%
disp('    Load sign restrictions')

% Unpack
N        = dims.N;
srScript = dataOpt.srScript;
shockVar0 = dataOpt.shockVar;
SRnew    = SR;

% patch for shockVar in char
shockVar = shockVar0;
if ischar(shockVar0)
shockVar = cellstr(shockVar0);
end

% Check that sign restrictions matrix is set up correctly
[c1,c2,c3] = size(SR);
if     c1 ~= N
    error(['ERROR: specify restrictions to all variables. Check ',srScript,'.m'])
elseif c2 ~= 3
    error(['ERROR: matrix SR set up incorrecty. Check ',srScript,'.m']')
elseif c3 ~= N
    
    if numel(shockVar) == 0
       error('ERROR: no shockVar defined. Specify as many shockVars as there are shocks you want to identify.') 
    end
    
    if c3 > N
    error(['ERROR: more shocks than endogenous variables. Check ',srScript,'.m']')
    elseif numel(shockVar) > c3
    error(['ERROR: more shocks in shockVar than what is defined in SR. Check ',srScript,'.m']')       
    end

end

if numel(shockVar) < c3
error(['ERROR: less shocks in shockVar than what is defined in SR. Check ',srScript,'.m']')    
end

dataOpt.nsr = c3;