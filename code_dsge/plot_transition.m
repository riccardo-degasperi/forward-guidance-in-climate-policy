%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Guidance in Climate Policy
% Authors: Riccardo Degasperi, Tara Hamadi, Filippo Natoli, Valerio Nispi Landi, Kevin Pallara
% % This file plots the IRF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clc; clear all; close all;


T=11;
T2=31;
SAVE=1;
load irf1            
load irf2
load irf3            

load GT
Tfull=length(y1)-1;         

Ta=1;
Tb=T;

% Define objects for plotting
FigW=26;FigH=13;
usedfont='Courier New';

cellvars=     {'pi',      'y',       'c',         'i',        'h',          'r',                 'e',        'price'};
cellvarlabels={'Inflation','Output', 'Consumption','Investment','Labor','Nominal rate',  'Emissions','Carbon price'};


           
%% All figures

figure;
set(gcf,'color', 'white',...
        'PaperUnits','centimeters','PaperSize',[6.5*4 6.5*2],...
        'PaperPosition',[0.25*4 0.25*2 5.75*4 5.75*2]);
set(gca,'FontName',usedfont);
tiledlayout(2,4)
j=1;
for j=1:8
nexttile
plottanda1=eval(['100*' cellvars{j} '1' ';']); 
plot([0:1:T-1],plottanda1(1:T),'b','Linewidth',1.5);
hold on
plottanda2=eval(['100*' cellvars{j} '2' ';']); 
plot([0:1:T-1],plottanda2(1:T),'k--','Linewidth',1.5);
hold on
plottanda4=eval(['100*' cellvars{j} '3' ';']); 
plot([0:1:T-1],plottanda4(1:T),'r:','Linewidth',1.5);
title(cellvarlabels{j},'Fontsize',11,'FontName',usedfont,'FontWeight','normal')
axis tight
ax = gca;
ax.FontName = usedfont;
if j==1 || j==6
ylabel('% points','Fontsize',10,'FontName',usedfont,'FontWeight','normal')
elseif j==8
ylabel('EUR per ton of Co2','Fontsize',10,'FontName',usedfont,'FontWeight','normal')
else
ylabel('%','Fontsize',10,'FontName',usedfont,'FontWeight','normal')
end
if j==8
xlabel('Horizon (quarters)','Fontsize',10,'FontName',usedfont,'FontWeight','normal')
end
j=j+1;   
end
legNames={'Baseline','Stance shock','Path shock'};    
hL = legend(legNames,'Fontsize',10,'FontName',usedfont,'FontWeight','normal');
hL.Orientation = 'horizontal';
hL.Layout.Tile = 'south';

if SAVE==1
if SPEC==0
print('-dpdf','-r500',['./figures/dsgebase.pdf'])
elseif SPEC==1
print('-dpdf','-r500',['./figures/dsgefast.pdf'])
elseif SPEC==2
print('-dpdf','-r500',['./figures/dsgehighprice.pdf'])
elseif SPEC==3
print('-dpdf','-r500',['./figures/dsgelowmu.pdf'])
elseif SPEC==4
print('-dpdf','-r500',['./figures/dsgeshort.pdf'])
end    
end

%% Carbon tax

figure;
set(gcf,'color', 'white',...
        'PaperUnits','centimeters','PaperSize',[13 13],...
        'PaperPosition',[0.25*2 0.25*2 5.75*2 5.75*2]);
set(gca,'FontName',usedfont);
j=8;
tiledlayout(1,1)
nexttile
plottanda1=eval(['100*' cellvars{j} '1' ';']); 
plot([0:1:T2-1],plottanda1(1:T2),'b','Linewidth',1.5);
hold on
plottanda2=eval(['100*' cellvars{j} '2' ';']); 
plot([0:1:T2-1],plottanda2(1:T2),'k--','Linewidth',1.5);
hold on
plottanda4=eval(['100*' cellvars{j} '3' ';']); 
plot([0:1:T2-1],plottanda4(1:T2),'r:','Linewidth',1.5);

title(cellvarlabels{j},'Fontsize',11,'FontName',usedfont,'FontWeight','normal')
axis tight
ax = gca;
ax.FontName = usedfont;
ylabel('EUR per ton of Co2','Fontsize',10,'FontName',usedfont,'FontWeight','normal')
xlabel('Horizon (quarters)','Fontsize',10,'FontName',usedfont,'FontWeight','normal')
legNames={'Baseline','Stance shock','Path shock'};    
hL = legend(legNames,'Fontsize',10,'FontName',usedfont,'FontWeight','normal');
hL.Orientation = 'horizontal';
hL.Layout.Tile = 'south';
if SAVE==1 && SPEC==0
print('-dpdf','-r500',['./figures/carbontax.pdf'])
end


