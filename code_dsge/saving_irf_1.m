%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Guidance in Climate Policy
% Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara
% This file saves the IRF in the benchmark scenario
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


dynare transition_shock

Tfull=length(y)-1;       
T=Tfull;


y1=y(1:T)/y(1)-1;
c1=(c(1:T))/(c(1))-1;
i1=(I(1:T))/(I(1))-1;
h1=(h(1:T))/(h(1))-1;
w1=(w(1:T))/(w(1))-1;
k1=(k(1:T))/(k(1))-1;
pi1=4*((pi(1:T))-(pi(1)));
r1=4*((r(1:T))-(r(1)));
piG1=4*((piG(1:T))-(piG(1)));
rreal1 = 4* (rreal(1:T) - rreal(1));
piB1=4*((piB(1:T))-(piB(1)));
X1=(X(1:T))/(X(1))-1;
e1=(e(1:T)/e(1))-1;
mu1=mu(1:T)-mu(1);
tau1=tt-tt(1);
price1=(price(1:T)-price(1))/100;
ab1=(nuM/(1+chi).*mu.^(1+chi))-(nuM/(1+chi).*mu(1).^(1+chi));
kG1=(kG(1:T))/(kG(1))-1;
kB1=(kB(1:T))/(kB(1))-1;

cW1=c;
hW1=h;

save irf1 y1 c1 i1 h1 k1 pi1 r1 piB1 X1 piG1 e1 mu1 tau1 ab1 price1 w1 rreal1 betta phi cW1 hW1 kG1 kB1
    %c_br1 rk_br1 w_br1 h_br1 y_br1 k_br1 q_br1 i_br1 r_br1 mc_br1 pi_br1 mu_br1 tau_br1 
