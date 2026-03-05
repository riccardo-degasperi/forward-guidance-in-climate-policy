%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Guidance in Climate Policy
% Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara
% This file computes the initial steady state
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear all;  
close all;
clc;
load par_trans

GTS=     [0 1 1];  % 0=benchmark, 1=the green transition is shocked (either path or stance)                           
SCENS=   [1 2 3];  % 0=baseline scenario, 2: stance shock, 3=path shock              
PLOTS=   [0 0 1];  % 0=dont' plot, 1=plot      
SPEC=0;            % 0=Figure 9; 1=D1; 2=D2; 3=D3; 4=D4; 


for kk=1:3
%% Setting
T=121;            % first period is the SS. Then in T-1 periods emissions go to 0
V=301;            % additional periods of simulation after T
Z=5;              % when the path/stance shock arrives (5 means that the shock arrives after 4 quarters, as the first period is the impact)
FAST=0;           % 0: green transition is linear, 1 green transition is convex (as in FeNL)
shock_AR=0.02;    % stance shock (200 basis points)
rho_AR=0.95;      % AR(1) parameter of stance shock
Zend=60;          % from Zend onward the path shock comes back to the baseline
shock_gradual=0.001; % how much the slope increases on impact after the path shock
step_gradual=0.001;  % how much the slope decreases after the shock in the convex transition


mu=1;             % final abatement
SCEN=SCENS(kk);    % 1: run saving_path1: benchmark rational
gtshock=GTS(kk);                
PLOT=PLOTS(kk);   % 1: run plot_transition; 0: no plot


clear y

%% Steady State
options = optimoptions('fsolve','MaxFunEvals',300000,'MaxIter',30000,'TolFun',1e-15);

rk=theta/betta-(1-delta);    
tau=nuM/nuE*mu^(chi);
r=piss*theta/betta;

load x0
x1=x0(1:2);

[x,u,Fval] = fsolve(@(xx) find_steady_clean(xx,alfa,phi,theta,delta,csi,g,zeta,chi,nuM,betta,nuE,tau,mu,varsig,epsG,epsB,rk),x1,options);
y=x(1);
pB=x(2);

if csi==1
pG=((1-zeta)^(1-zeta)*zeta^(zeta)*(pB)^(-zeta))^(1/(1-zeta));
else
pG=(1/(1-zeta)*(1-zeta*(pB)^(1-csi)))^(1/(1-csi));
end
mcG=pG*(epsG-1)/epsG;            
mcB=pB*(epsB-1)/epsB-tau*(1-mu)*nuE-nuM/(1+chi)*mu^(1+chi); 
yB=zeta*((pB)^(-csi)*y);
yG=(1-zeta)*(pG)^(-csi)*y;   
kG=alfa*yG*mcG/(rk);
kB=alfa*yB*mcB/(rk);                       
hB=(yB/((kB)^(alfa)))^(1/(1-alfa));
hG=(yG/((kG)^(alfa)))^(1/(1-alfa));
k=theta*(kG+kB);
i=(1-(1-delta)/theta)*k;                           
w=(1-alfa)*yG*mcG/hG;
h=hB+hG;
c=y-i-g-nuM/(1+chi)*mu^(1+chi)*yB;  
lam=(theta-betta*varsig)/(c*(theta-varsig));
e=(1-mu)*nuE*yB;
tau_clean=tau;
y_clean=y;
mu_clean=mu;
e_clean=e;
pB_clean=pB;




%% Build tax process
% Benchmark
if FAST==0 % linear transition
step=(tau_clean-tau_start)/(T-1); 
tt=zeros(T+V+1,1);
tt(1)=tau_start;
for j=2:1:T
    tt(j)=tt(j-1)+step;
end
tt(T+1:T+V+1)=ones(V+1,1)*tau_clean;

else % convex transition
S=tau_clean-tau_start;
step_step=2*(S)/((T-1)*(T));
step=zeros(T,1);
step(T)=0;
for j=2:1:T
    step(T-j+1)=step(T-j+2)+step_step;
end
tt=zeros(T+V,1);
tt(1)=tau_start;
for j=2:1:T+1
    tt(j)=tt(j-1)+step(j-1);
end
tt(T+2:T+1+V)=ones(V,1)*tau_clean;
end


% Stance shock
tt_shock_ar=tt;

shock_ar=zeros(length(tt)-(Z-1),1);
shock_ar(1)=shock_AR;
for jj=2:length(shock_ar)
shock_ar(jj)=rho_AR*shock_ar(jj-1);
end
tt_shock_ar(Z:end)=tt(Z:end)+shock_ar;

% Ensure that the tax is never higher than the maximum level (the level such that emissions are 0)
for j=1:length(tt_shock_ar)
 if tt_shock_ar(j)>tau_clean
     tt_shock_ar(j)=tau_clean;
 end
end



% Path shock
tt_slope = tt;
original_slope = (tau_clean - tau_start) / (T - 1);  % original slope

if FAST==0

tt_slope(Z) = tt(Z-1) + original_slope + shock_gradual;

for j = Z+1:Zend
    scale_factor = 1 - (j - Z) / (Zend - Z);    % how much the slope decreases
    new_slope = original_slope + shock_gradual * scale_factor;  % new slope after the shock
    tt_slope(j) = tt_slope(j-1) + new_slope;
end

% After Zend periods, the tax comes back to the original shock
for j = Zend+1:length(tt)
    if j <= T
        % Back to the original shock
        tt_slope(j) = tt_slope(j-1) + original_slope;
    else
        % After the last period of the transition, the tax is constant
        tt_slope(j) = tt_slope(j-1);
    end
end

else    

tt_slope(Z) = tt(Z) + shock_gradual;

for j = Z+1:Zend
    scale_factor = 1 - (j - Z) / (Zend - Z);
    if j-1 <= T
        tt_slope(j) = tt_slope(j-1) + step(j-1) + step_gradual * scale_factor;
    else
        tt_slope(j) = tt_slope(j-1) + step_gradual * scale_factor; 
    end
end

for j = Zend+1:length(tt)
    if j-1 <= T
        tt_slope(j) = tt_slope(j-1) + step(j-1);
    else
        tt_slope(j) = tt_slope(j-1); 
    end
end
end

% Ensure that the tax is never higher than the upper bound
tt_slope(tt_slope > tt(end)) = tt(end);


%% Save parameters
save par_clean tau_clean y_clean tt mu_clean e_clean pB_clean gtshock tt_shock_ar
save GT gtshock Z SPEC


if SCEN==1
tt_shock=NaN;
saving_irf_1;
elseif SCEN==2
tt_shock=tt_shock_ar;
saving_irf_2
elseif SCEN==3
tt_shock=tt_slope;
saving_irf_3
end
if PLOT==1
    plot_transition
end

end