%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Guidance in Climate Policy
% Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara
% This file saves the IRF under the stance shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



dynare transition_shock

Tfull=length(y)-1;         % path length
T=Tfull;
c=yy(strmatch('c',M_.endo_names,'exact'),:);
w =yy(strmatch('w',M_.endo_names,'exact'),:);
h =yy(strmatch('h',M_.endo_names,'exact'),:);
y =yy(strmatch('y',M_.endo_names,'exact'),:);
k =yy(strmatch('k',M_.endo_names,'exact'),:); 
q =yy(strmatch('q',M_.endo_names,'exact'),:);
I =yy(strmatch('I',M_.endo_names,'exact'),:);
r =yy(strmatch('r',M_.endo_names,'exact'),:);
pi=yy(strmatch('pi',M_.endo_names,'exact'),:);
mu=yy(strmatch('mu',M_.endo_names,'exact'),:);
e =yy(strmatch('e',M_.endo_names,'exact'),:);
X =yy(strmatch('X',M_.endo_names,'exact'),:);
piG=yy(strmatch('piG',M_.endo_names,'exact'),:);
piB=yy(strmatch('piB',M_.endo_names,'exact'),:);
price=yy(strmatch('price',M_.endo_names,'exact'),:);
rreal=yy(strmatch('rreal',M_.endo_names,'exact'),:);
kG=yy(strmatch('kG',M_.endo_names,'exact'),:);
kB=yy(strmatch('kB',M_.endo_names,'exact'),:);


y2=y(1:T)/y(1)-1;
c2=(c(1:T))/(c(1))-1;
i2=(I(1:T))/(I(1))-1;
h2=(h(1:T))/(h(1))-1;
w2=(w(1:T))/(w(1))-1;
k2=(k(1:T))/(k(1))-1;
pi2=4*((pi(1:T))-(pi(1)));
r2=4*((r(1:T))-(r(1)));
piG2=4*((piG(1:T))-(piG(1)));
rreal2 = 4* (rreal(1:T) - rreal(1));
q2=q(1:T)/q(1)-1;
X2=(X(1:T))/(X(1))-1;
e2=(e(1:T)/e(1))-1;
mu2=mu(1:T)-mu(1);
tau2=tt-tt(1);
price2=(price(1:T)-price(1))/100;
ab2=(nuM/(1+chi).*mu.^(1+chi))-(nuM/(1+chi).*mu(1).^(1+chi));
piB2=4*((piB(1:T))-(piB(1)));
cW2=c;
hW2=h;
kG2=(kG(1:T))/(kG(1))-1;
kB2=(kB(1:T))/(kB(1))-1;

save irf2 y2 c2 i2 h2 k2 pi2 r2 q2 X2 piG2 e2 mu2 tau2 ab2 price2 w2 rreal2 cW2 hW2 piB2 kG2 kB2
    