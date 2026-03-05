%-------------------------------------------------------------------------%
% MATLAB TOOLBOX                                                          %
% By Riccardo Degasperi                                                   %
%                                                                         %
% riccardo.degasperi@bancaditalia.it                                      %
%-------------------------------------------------------------------------%

addpath([pwd '/subroutines/'])                   %path to child functions
addpath([pwd '/data/'])                          %path to data folder

%-------------------------------------------------------------------------%
% Data import
disp('--> Preliminary operations...')

[data,dataOpt] = dataLoaderBVAR(modelOpt,dataOpt,miscOpt);

if plotOpt.exploreSeries

    lags = 6;   % select the number of lags for ADF
    seriesExplorer(data,dataOpt,plotOpt,lags)

end

%-------------------------------------------------------------------------%
% Checks and patches
checks


%-------------------------------------------------------------------------%
% Plot series
if plotOpt.plotTimeSeries

    plotSeries(data,dataOpt,plotOpt)
    
end


%-------------------------------------------------------------------------%
% Make X and Y matrices
[data,dims,hp,dataOpt] = makeXY(data,hp,modelOpt,dataOpt);
plotOpt.ylab = dataOpt.ylab;


%-------------------------------------------------------------------------%
% Load sign and magnitude restrictions
if strcmp(modelOpt.identification,'sign') || strcmp(modelOpt.identification,'jarocinski-karadi') || strcmp(modelOpt.identification,'iv+sign')
    
    eval(dataOpt.srScript);
    modelOpt.SR = SRaugment(modelOpt.SR,dims,dataOpt);

    if ~isempty(dataOpt.mrScript)
        eval(dataOpt.mrScript);
        modelOpt.MR = MRaugment(modelOpt.MR,dims,dataOpt);
    else
        modelOpt.MR = zeros(size(modelOpt.SR,1),5,size(modelOpt.SR,3));
    end
    
end



%-------------------------------------------------------------------------%
% Full Gibbs sampler for Bayesian VAR (dummies or triangular)
if modelOpt.fullGibbs

out = fullGibbsSampler(data,dims,hp,modelOpt,dataOpt,miscOpt);

% Plot IRFs
chartName = [char(dataOpt.units),'_FullGibbs'];
BVARplots(out.IRFsulm,chartName,dataOpt,plotOpt);


%-------------------------------------------------------------------------%
% Full Gibbs sampler for in-sample structural scenario analysis
% (Antolin-Diaz, Petrella, Rubio-Ramirez, 2021; Breitenlechner, Georgiadis, Schumann, 2022)
elseif modelOpt.SSAfullIRF

    SSAout = SSAfullGibbsIRF(data,dims,hp,modelOpt,dataOpt,miscOpt,plotOpt,ssaOpt);

else
%-------------------------------------------------------------------------%
% Estimating BVAR via dummy observations
% (Banbura, Giannone, Reichlin, 2010; Giannone, Lenza, Primiceri, 2015)
if modelOpt.dummies
    
    % Estimate BVAR
    BVAR_d = BVARestimate_dummies(data,dims,hp,modelOpt,dataOpt,miscOpt,plotOpt);
    
    % Structural identification
    BVAR_d = BVARidentification(BVAR_d,data,dims,modelOpt,dataOpt,miscOpt);
    
    % Plot IRFs
    chartName = [char(dataOpt.units),'_DummyObs'];
    BVARplots(BVAR_d.IRFsulm,chartName,dataOpt,plotOpt);

    % Structural Scenario Analysis (IRF)
    if miscOpt.SSAirf
        plotOpt.customLabels = {};
        plotOpt.estimation = 'DummyObs';
        BVAR_d.SSA = getSSA_irf(BVAR_d,data,dims,dataOpt,plotOpt,ssaOpt);
    end
    
end


end


