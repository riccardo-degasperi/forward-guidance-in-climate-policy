%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Guidance in Climate Policy
% Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara
% This file simulates the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all;
warning off



load GT
  
%%%%%%%%%%%%%%%%%%%%%%%Endogenous Variables %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var
c
w 
h
y 
k
q
I
r
pi
mu
e
X
piG
piB
price
rreal
lambda
rk
kG
kB
hG
hB
yG
yB
mcG
mcB
pG
pB
tau_fake

;
%%
%%%%%%%%%%%%%%%%%%%%%%%Exogenous Variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varexo 
tau
;  
    
%%
%%%%%%%%%%%%%%%%%%%%%%%Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

parameters

betta 
alfa 
epsG
epsB
delta
phi
chi
theta
deltax
nuM
nuE
y_start
mu_start
X_start
rss
piss
kappaI
phipi
phiy
rhor
kappaPG
kappaPB
zeta
varsig
csi
S1
S2
S3
tau_start
mu_clean
y_clean
pB_start
e_start
erow
iota
g
;

load par_trans;  % load mat file created in console_ss
load par_clean;  % load mat file created in console_ss_clean

for jj=1:length(M_.param_names)
set_param_value(M_.param_names{jj},eval(M_.param_names{jj})); 
end;


%%
%%%%%%%%%%%%%%%%%%%%%%%Non-Linear Model%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model;

%Households
h^(phi)=w*lambda;                                                      %(1)                                                         
lambda=1/(c-varsig*c(-1)/theta)-betta*varsig/(c(1)*theta-varsig*c);    %(2)
1=betta*lambda(1)/(lambda*theta)*r/pi(1);                              %(3)                   
1=betta*lambda(1)/(lambda*theta)*(rk(1)+(1-delta)*q(1))/q;             %(4)                
1=q*(1-kappaI/2*(I/I(-1)*theta-theta)^2-kappaI*(I/I(-1)*theta-theta)*I/I(-1)*theta)+kappaI*betta*lambda(1)/(lambda*theta)*q(1)*(I(1)/I*theta-theta)*(I(1)/I*theta)^2; %(5)   
k=(1-delta)*k(-1)/theta+(1-kappaI/2*(I/I(-1)*theta-theta)^2)*I;        %(6)
% Final good firm
y^((csi-1)/csi)=(1-zeta)^(1/csi)*(yG)^((csi-1)/csi)+(zeta)^(1/csi)*(yB)^((csi-1)/csi);  %(7)
yG=(1-zeta)*(pG)^(-csi)*y;                                             %(8)
yB=zeta    *(pB)^(-csi)*y;                                             %(9)
% Intermediated Firms
(piG-piG(-1)^(iota)*piss^(1-iota))*piG=betta*lambda(1)/lambda*pG(1)*yG(1)/(pG*yG)*piG(1)*(piG(1)-piG^(iota)*piss^(1-iota))+epsG/kappaPG*(mcG/pG-(epsG-1)/epsG); %(10)
pG/pG(-1)=piG/pi;                                                     %(11)
yG=(kG)^(alfa)*hG^(1-alfa);                                           %(12)
(1-alfa)*mcG*yG=w*hG;                                                 %(13)
alfa*mcG*yG=rk*kG;                                                    %(14)
(piB-piB(-1)^(iota)*piss^(1-iota))*piB=betta*lambda(1)/lambda*pB(1)*yB(1)/(pB*yB)*piB(1)*(piB(1)-piB^(iota)*piss^(1-iota))+epsB/kappaPB*((mcB+tau*(1-mu)*nuE+nuM/(1+chi)*mu^(1+chi))/pB-(epsB-1)/epsB); %(15)
pB/pB(-1)=piB/pi;                                                     %(16)
yB=(kB)^(alfa)*hB^(1-alfa);                                           %(17)
(1-alfa)*mcB*yB=w*hB;                                                 %(18)
alfa*mcB*yB=rk*kB;                                                    %(19)
mu=(nuE/nuM*tau)^(1/chi);                                             %(20)
% Pollution
X=(1-deltax)*X(-1)/theta+e+erow;                                      %(21)
e=(1-mu)*nuE*yB;                                                      %(22)
% Market clearing 
k(-1)/theta=kG+kB;                                                    %(23)
h=hG+hB;                                                              %(24)
y=c+I+g+yB*nuM/(1+chi)*mu^(1+chi)+pB*yB*kappaPB/2*(piB-piB(-1)^(iota)*piss^(1-iota))^2+pG*yG*kappaPG/2*(piG-piG(-1)^(iota)*piss^(1-iota))^2;   %(25)
% Policy
r/rss=((pi/piss)^(phipi)*(y/y(-1))^(phiy))^(1-rhor)*(r(-1)/rss)^(rhor);    %(26)
% Other defininitions
rreal = r / pi(1);                                                    %(27)                                                              
price=tau*S1*S2/S3;                                                   %(28)
tau_fake=tau;                                                         %(29)
end;

initval;
pi=piss;
piG=pi;
piB=pi;
q=1;
r=pi*theta/betta;             
rk=theta/betta-(1-delta); 
y=y_start;
pB=pB_start;
e=e_start;
mu=mu_start;
tau=(nuM/nuE)*(mu^(chi));
X=X_start;
pG=(1/(1-zeta)*(1-zeta*(pB)^(1-csi)))^(1/(1-csi));
mcG=pG*(epsG-1)/epsG;            
mcB=pB*(epsB-1)/epsB-tau*(1-mu)*nuE-nuM/(1+chi)*mu^(1+chi); 
yB=zeta*((pB)^(-csi)*y);
yG=(1-zeta)*(pG)^(-csi)*y;   
kG=alfa*yG*mcG/(rk);
kB=alfa*yB*mcB/(rk);                       
hB=(yB/((kB)^(alfa)))^(1/(1-alfa));
hG=(yG/((kG)^(alfa)))^(1/(1-alfa));
k=theta*(kG+kB);
I=(1-(1-delta)/theta)*k;                           
w=(1-alfa)*yG*mcG/hG;
h=hB+hG;
c=y-I-g-nuM/(1+chi)*mu^(1+chi)*yB;  
lambda=(theta-betta*varsig)/(c*(theta-varsig));
price=tau*S1*S2/S3;
rreal = r/pi;
tau_fake=tau;
end;

steady;

endval;
pi=piss;
piG=pi;
piB=pi;
q=1;
r=pi*theta/betta;             
rk=theta/betta-(1-delta); 
y=y_clean;
pB=pB_clean;
tau=tau_clean;
mu=(nuE/nuM*tau)^(1/chi);  
pG=(1/(1-zeta)*(1-zeta*(pB)^(1-csi)))^(1/(1-csi));
mcG=pG*(epsG-1)/epsG;            
mcB=pB*(epsB-1)/epsB-tau*(1-mu)*nuE-nuM/(1+chi)*mu^(1+chi);
yB=zeta*((pB)^(-csi)*y);
yG=(1-zeta)*(pG)^(-csi)*y;   
kG=alfa*yG*mcG/(rk);
kB=alfa*yB*mcB/(rk);                       
hB=(yB/((kB)^(alfa)))^(1/(1-alfa));
hG=(yG/((kG)^(alfa)))^(1/(1-alfa));
k=theta*(kG+kB);
I=(1-(1-delta)/theta)*k;                           
w=(1-alfa)*yG*mcG/hG;
h=hB+hG;
c=y-I-g-nuM/(1+chi)*mu^(1+chi)*yB;  
lambda=(theta-betta*varsig)/(c*(theta-varsig));
price=tau*S1*S2/S3;
rreal = r/pi;
tau_fake=tau;
e=(1-mu)*nuE*yB;
X=(e+erow)/(1-(1-deltax)/theta);

end;
 
steady;

t2=tt(2:end);    
shocks;
var tau;
periods 1:422;
values (t2);
end;


perfect_foresight_setup(periods=422);
perfect_foresight_solver; 
yy=oo_.endo_simul(:,1:Z-1);     % We take only the first Z columns of simulations: the steady state and the Z-1 initial periods

% if green transition shock
if gtshock==1
oo_.exo_simul(1:end-(Z-1))=tt_shock(Z-1:end);  % In period Z the green transition changes
oo_.endo_simul(:,1)=yy(:,end);
perfect_foresight_solver; 
yy=[yy,oo_.endo_simul(:,2:end)];
end