%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Guidance in Climate Policy
% Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara
% This function computes the final steady state. It is called in console_ss_clean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ F ]  = find_steady_clean(x,alfa,phi,theta,delta,csi,g,zeta,chi,nuM,betta,nuE,tau,mu,varsig,epsG,epsB,rk)
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

F=[w*hB-(1-alfa)*yB*mcB
    h^(phi)-w*lam
    ];

end

