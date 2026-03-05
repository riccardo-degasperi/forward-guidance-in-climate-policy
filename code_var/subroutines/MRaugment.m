function MRnew = MRaugment(MR,dims,dataOpt)
% MRaugment checks that the magnitude restriction script has been compiled
% properly and prepares the restrictions for later use.
%
% last modified: 26/03/2024
% riccardo.degasperi@bancaditalia.it
%-------------------------------------------------------------------------%
disp('    Load magnitude restrictions')

% Unpack
N        = dims.N;
mrScript = dataOpt.mrScript;
shockVar0 = dataOpt.shockVar;
MRnew    = MR;

% patch for shockVar in char
shockVar = shockVar0;
if ischar(shockVar0)
shockVar = cellstr(shockVar0);
end

% Check that sign restrictions matrix is set up correctly
[c1,c2,c3] = size(MR);
if     c1 ~= N
    error(['ERROR: specify restrictions to all variables. Check ',mrScript,'.m'])
elseif c2 ~= 5
    error(['ERROR: matrix MR set up incorrecty. Check ',mrScript,'.m']')
elseif c3 ~= N

    if numel(shockVar) == 0
       error('ERROR: no shockVar defined. Specify as many shockVars as there are shocks you want to identify.') 
    end
    
    if c3 > N
    error(['ERROR: more shocks than endogenous variables. Check ',mrScript,'.m']')
    elseif numel(shockVar) > c3
    error(['ERROR: more shocks in shockVar than what is defined in MR. Check ',mrScript,'.m']')    
    end

end

if numel(shockVar) < c3
error(['ERROR: less shocks in shockVar than what is defined in MR. Check ',mrScript,'.m']')    
end

dataOpt.nmr = c3;