%-------------------------------------------------------------------------%
% F O R W A R D   G U I D A N C E   I N   C L I M A T E   P O L I C Y     %
%                                                                         %
% By Riccardo Degasperi, Tara Hamadi, Filippo Natoli,                     %
% Valerio Nispi Landi and Kevin Pallara                                   %
%                                                                         %
% Replica Figures C.1 and C.2                                             %
%-------------------------------------------------------------------------%
clear all
close all
clc

addpath(genpath(pwd))  %path to child functions

proxies = {'F1_EUA';'F2_EUA'};
UU = size(proxies,1);

for z = 1:UU

clearvars -except z UU proxies
close all

%-------------------------------------------------------------------------%
% Type of model
modelOpt.fullGibbs      = 0;
modelOpt.dummies        = 1;                    %1:on; 0:off

% Model details
modelOpt.p              = 6;                   %number of lags for endogenous variables
modelOpt.pexo           = 0;                    %lags for exogenous variables (must be < p)
modelOpt.h              = 36;                   %horizons
modelOpt.constant       = 1;                    %0:no constant; 1:constant; 2:constant+trend
modelOpt.sim            = 1600;                 %number of simulations
modelOpt.burnin         = 600;                  %burnin
modelOpt.jump           = 1;                    %stores a draw every jump draws

% Date cutoffs
dataOpt.lT              = 'full';               %lower date cutoff (either 'full' or a date in text format)
dataOpt.uT              = 'full';               %upper date cutoff (either 'full' or a date in text format)
dataOpt.frequency       = 'monthly';            %frequency of data, can be 'daily', 'weekly', 'monthly', 'quarterly', 'yearly' (also displayed in the charts)


% Series selection
dataOpt.dataset         = 'data_replica.xlsx';  %dataset name including extension
dataOpt.shockVar        = {'EKESCPENF_SA'};     %shock of interest (position in the impact matrix)
dataOpt.shockUnit       = {};                   %unit originating the shock (leave empty if policy variable is in sheet "Global")
dataOpt.shockSign       = 'positive';           %either 'positive' or 'negative'
dataOpt.shockSize       = '1sd';                %magnitude of the shocks. Either '1sd' or a scalar (also negative)
dataOpt.varendo         = {'EKESCPENF_SA','GHGTotal_IP','EKCPHARMF_SA','EKIPTOTG','EMINTER3','EKESUNEMO','DJSTO50','ROILEU'};
dataOpt.varendo_global  = {};
dataOpt.varplot         = dataOpt.varendo;
dataOpt.varexo          = {};                   %exogenous variables
dataOpt.units           = {'data_kaenzig'};       %units (can be more than one)

% Identification
modelOpt.identification = 'iv';                 %'cholesky', 'iv', 'sign', 'jarocinski-karadi', 'iv+sign'
modelOpt.inference      = '';                   %for frequentist only: 'montecarlo','bootstrap','wild bootstrap'
modelOpt.maxSRtries     = 1000;                 %maximum number of tries for each iteration of the sign restriction
dataOpt.srScript        = 'SRgreen';            %name of script defining sign restrictions (customise it before launch)
dataOpt.mrScript        = '';                   %name of script defining magnitude restrictions (customise it before launch)
dataOpt.variv           = proxies(z);
dataOpt.shockIV         = [1];                  %position of shocks in shockVar to be identified with proxies

% Options for BVAR
hp.GLP                  = 1;                    %optimal prior selection
hp.hyperpriors          = 1;                    %(if GLP=1): 1:uses hyperpriors; 0:ML only
hp.niw                  = 1;                    %Normal Inverse Wishart prior
hp.soc                  = 0;                    %Sum of Coefficients prior
hp.cop                  = 0;                    %Copersistence prior
hp.lag                  = 0;                    %optimise lag-decay prior
hp.std                  = 0;                    %optimise diagonal elements of the scale matrix of the IW prior on the residual variance
hp.cross                = 0;                    %optimise extra shrinkage parameter on cross-variable coefficients
hp.covid                = 0;                    %pandemic priors (Cascaldi-Garcia, 2024)

% Setting for pandemic priors time dummies
dataOpt.lT_covid        = '2020m3';             %first time dummy for pandemic priors
dataOpt.uT_covid        = '2020m8';             %last time dummy for pandemic priors

% Initial values for priors
hp.l1                   = 0.2;                  %Normal-Inverse-Wishart
hp.l2                   = 1;                    %lag decay (=1 is off)
hp.l3                   = 1;                    %Sum-of-coefficients (=1 is off)
hp.l4                   = 1;                    %Co-persistence (=1 is off)
hp.l5                   = 1000;                 %constant and exogenous
hp.l6                   = 1;                    %extra shrinkage on cross-variable parameters (=1 is off)
hp.l7                   = 0.05;                 %informativeness of prior for pandemic time dummies
hp.eps                  = 1e-10;                %dogmatic shrinkage parameter

% Options for Structural Scenario Analysis
modelOpt.SSAfullIRF     = 0;
modelOpt.SSAfullTS      = 0;
miscOpt.SSAirf          = 0;                    %1:on; 2:off
miscOpt.SSAts           = 0;                    %1:on; 2:off
ssaOpt.SSAlT            = '';                   %(for SSAts) sets lower Xlim for scenarios plot
ssaOpt.SSAscript        = '';                   %name of function defining SSA restrictions (customise it before launch)
ssaOpt.KLtresh          = 0;                    %threshold for KL divergence
ssaOpt.impulse          = {};                   %generate IRFs to shocks in position 'impulse' (you can set this to dataOpt.shockVar)
ssaOpt.SSAy             = {};                   %variables whose path is restricted (can define many in this format SSAy = {{'var1','var2'};{'var1'};{'var1','var2','var3'}};)
ssaOpt.SSAe             = {};                   %variables whose innovations are used to deliver the path
ssaOpt.customLabels     = {};

% Charts options
plotOpt.saveCharts      = 1;                    %1:save graphs to folder ./plots
plotOpt.saveFig         = 0;
plotOpt.exploreSeries   = 0;                    %one-by-one plots, ACF, PACF and ADF of series in levels and first differences
plotOpt.plotTimeSeries  = 0;                    %1:plot time series; 0:don't
plotOpt.comparisons     = 0;                    %(if Triangular) plots comparisons Triangular vs. NIW
plotOpt.plotType        = 'singleShock';        %either 'singleShock' or 'allShocks'
plotOpt.whichPlot       = 'varendo';            %either 'all','varendo','varplot','colour'
plotOpt.cb90            = 1;                    %1:plots also 90% confidence bands; 0:plots only 68% confidence bands
plotOpt.folder          = ['figures_C1_C2/',dataOpt.variv{:}];  %either a folder name or []

plotOpt.allSeries.plotRows = 2;                 %max number of rows per plot
plotOpt.allSeries.plotCols = 4;
plotOpt.varplot.plotRows   = 4;                 %max number of rows per plot
plotOpt.varplot.plotCols   = 4;
plotOpt.font               = 'Courier New';     %font for charts

% Misc options
miscOpt.dropExplosive   = 1;                    %drop explosive draws
miscOpt.dateCheck       = 0;                    %check sample length when option 'full' is on. [Requires user input]
miscOpt.interpolate     = 1;                    %1:interpolate NaNs; 0:don't
miscOpt.interpolCase    = 'spline';             %'MA':centered MA; 'spline':cubic spline
miscOpt.ws              = 1;                    %window size (for 'MA' option)


%-------------------------------------------------------------------------%
% Run analysis
main

end