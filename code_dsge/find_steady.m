%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Guidance in Climate Policy
% Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara
% This function computes the initial steady state. It is called in console_ss
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [ F ]  = find_steady(x,phi,theta,delta,csi,G,zeta,chi,nuM,mu,rk,I_y,gdpEUR,XGtC,price,S3,RoW,deltax,betta,epsG,epsB,varsig)
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

F=[w*hB-(1-alfa)*yB*mcB
    h^(phi)-w*lam
    e-(1-mu)*nuE*yB
    ];

end

