%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Guidance in Climate Policy
% Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara
% This file saves the IRF of the path shock
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



dynare transition_shock

Tfull=length(y)-1;        
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

y3=y(1:T)/y(1)-1;
c3=(c(1:T))/(c(1))-1;
i3=(I(1:T))/(I(1))-1;
h3=(h(1:T))/(h(1))-1;
piG3 = 4*(piG(1:T)  - piG(1));
w3=(w(1:T))/(w(1))-1;
k3=(k(1:T))/(k(1))-1;
pi3=4*((pi(1:T))-(pi(1)));
r3=4*((r(1:T))-(r(1)));
rreal3 = 4* (rreal(1:T) - rreal(1));
q3=q(1:T)/q(1)-1;
X3=(X(1:T))/(X(1))-1;
e3=(e(1:T)/e(1))-1;
mu3=mu(1:T)-mu(1);
tau3=tt-tt(1);
price3=(price(1:T)-price(1))/100;
ab3=(nuM/(1+chi).*mu.^(1+chi))-(nuM/(1+chi).*mu(1).^(1+chi));
piB3=4*((piB(1:T))-(piB(1)));

save irf3 y3 c3 i3 h3 k3 pi3 r3 piG3 q3 X3 e3 mu3 tau3 ab3 price3 w3 rreal3 piB3
