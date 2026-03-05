%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Guidance in Climate Policy
% Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara
% This file computes the initial steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; 
clc; 
close all;

%% Important choices
phipi=2.74;              % mp response to inflation
phiy=0;                  % response to output growth
rhor=0.93;               % monetary policy inertia (0.86 previous calibration)
calvoG=0.82;             % Calvo parameters for green goods (NAWM)
fmG=0.1418;              % frequency of monthly price adjustment for G according to Del Negro et al.
fmB=0.205;               % frequency of monthly price adjustment for B according to Del Negro et al.
fqG=1-(1-fmG)^3;         % frequency of quarterly price adjustment for G according to Del Negro et al.
fqB=1-(1-fmB)^3;         % frequency of quarterly price adjustment for B according to Del Negro et al.
calvoB=1-fqB/fqG*(1-calvoG); % Calvo parameters for B, to get a flat PC but to respect the proportion between adjustment in B and G

%% Other parameters
epsG=3.8571;                  % elasticity of substitution btw differentiated G goods
epsB=3.8571;                  % elasticity of substitution btw differentiated B goods
theta=1.0038;                 % steady-state growth
delta=0.025;                  % depreciation rate firms
phi=2;                        % inverse of Frisch elasticity
csi=2;                        % elasticity of substitution btw wgreen and brown production (CHM 21)
pi=1.005;                     % steady-state inflation. 
deltax=0.0035;                % pollution depreciation
zeta=0.8;                     % brown sector size (GKRW 21)
chi=1.6;                      % curvature the abatement function
nuM=0.074*(1+chi);            % abatement shifter
gdpEUR=3022.4;                % euro area quarterly GDP in 2019Q4 (mld EUR) https://appsso.eurostat.ec.europa.eu/nui/submitViewTableAction.do
XGtC=870.1476;                % pollution in Gigatons of Carbon in 2018 (https://ourworldindata.org/atmospheric-concentrations, 1ppm=2.13 GtC)
varsig=0.8;                   % habits formation (to match volatility)

%% Targets
rr=1+0.02/4;                  % real interest rate
I_y=0.21;                     % investment/GDP
G=0.215;                      % public spending/GDP ratio
RoW=15.31;                    % RoW emission / country emission (see excel file)
mu=0;                         % fraction of emissions abated
price=65;                     % carbon price in Euro under full abatement  
S3=3.67;                      % 1 unit of carbon=3.67 units of CO2

%% Steady State

betta=theta/rr;                        % discount factor
r=rr*pi;                               %  nominal interest rate
rk=theta/betta-(1-delta);              % rental rate of capital

load x0
x0=x0+0.001;
options = optimoptions('fsolve','MaxFunEvals',300000,'MaxIter',30000,'TolFun',1e-15);

[x,F,Fval] = fsolve(@(xx) find_steady(xx,phi,theta,delta,csi,G,zeta,chi,nuM,mu,rk,I_y,gdpEUR,XGtC,price,S3,RoW,deltax,betta,epsG,epsB,varsig),x0,options);
y=x(1);
pB=x(2);
e=x(3);


if csi==1
pG=((1-zeta)^(1-zeta)*zeta^(zeta)*(pB)^(-zeta))^(1/(1-zeta));
else
pG=(1/(1-zeta)*(1-zeta*(pB)^(1-csi)))^(1/(1-csi));
end
erow=RoW*e;
X=(e+erow)/(1-(1-deltax)/theta); 
S1=gdpEUR/y;                         
S2=XGtC/X;  
nuE=nuM/price*S1*S2/S3;
tau=nuM/nuE*mu^(chi);
mcG=pG*(epsG-1)/epsG;            
mcB=pB*(epsB-1)/epsB-tau*(1-mu)*nuE-nuM/(1+chi)*mu^(1+chi); 
yB=zeta*((pB)^(-csi)*y);
yG=(1-zeta)*(pG)^(-csi)*y;  
i=I_y*y;
k=i/(1-(1-delta)/theta); 
alfa=k/(theta*(yG*mcG/rk+mcB*yB/rk));
kG=alfa*yG*mcG/(rk);
kB=alfa*yB*mcB/(rk);                       
hB=(yB/((kB)^(alfa)))^(1/(1-alfa));
hG=(yG/((kG)^(alfa)))^(1/(1-alfa));
w=(1-alfa)*yG*mcG/hG;
h=hB+hG;
g=G*y;
c=y-i-g-nuM/(1+chi)*mu^(1+chi)*yB;  
lam=(theta-betta*varsig)/(c*(theta-varsig));


if Fval>-1           
x0=[y,pB,e];
save x0 x0
end


%% Other parameters that do not affect the ss

kappaI=10.78;              % investment adjustment cost (as in CEE). If 0, q is constant
kappaPG=(epsG-1)*calvoG/(pi^2*(1-calvoG)*(1-betta*calvoG)); 
kappaPB=(epsB-1)*calvoB/(pi^2*(1-calvoB)*(1-betta*calvoB)); 
iota=0.0;                 % inflation inertia
X=(e+erow)/(1-(1-deltax)/theta);

%% Steady state values
piss=pi;
y_start=y;
mu_start=mu;
X_start=X;
tau_start=tau;
e_start=e;
pB_start=pB;
rss=r;


%% Save parameters

save par_trans betta alfa  delta phi chi theta deltax erow g nuE nuM...
               rss piss y_start mu_start X_start tau_start e_start phiy...
               kappaI phipi rhor kappaPG S1 S2 S3 tau_start iota G g...
               varsig zeta csi pB_start epsG epsB rss calvoB kappaPB calvoG kappaPG

console_ss_clean
