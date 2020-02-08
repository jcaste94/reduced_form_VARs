%==========================================================================
%           ECONOMETRICS IV - SPRING 2020 - Prof. Schorfehide 
%                       REDUCED-FORM VARs                     
%
%
% Author: Juan Castellanos Silvan (from Luigi Boccola's original code)
% Date  : 02/06/2020
%==========================================================================


%=========================================================================
%                             HOUSEKEEPING
%=========================================================================

close all
clear all
clc

%=========================================================================
%       DIRECT MONTE CARLO SAMPLING FROM POSTERIOR OF VAR PARAMETERS
%=========================================================================

load("resultsVAR")

%=========================================================================
%           FIGURE 1: EVOLUTION OF ACTUAL DATA
%=========================================================================

% Attributes 
linesty={'-','-','-','-'};
color={'k','b','r','g'};
marker={'none','none','none','none'};
legend_list = {'Output', 'Inflation', 'Interest rate', 'Inverse velocity'};

% xaxis
time = 1966:0.25:2006;

figure(1);
p = plot(time, YYact, 'linewidth', 1);
size_p=size(p);
grid on;
for i=1:size_p(1);
    set(p(i),'LineStyle',linesty{i},'Color',color{i},'Marker',marker{i});
    xlabel("Time")
end
h=legend(legend_list, 'orientation', 'horizontal');
set(h,'Fontsize',10);
set(h,'Position',[0.45 0.01 0.1 0.05]);

x = 29.7;                  % A4 paper size
y = 21.0;                  % A4 paper size
xMargin = 1;               % left/right margins from page borders
yMargin = 1;               % bottom/top margins from page borders
xSize = x - 2*xMargin;     % figure size on paper (widht & hieght)
ySize = y - 2*yMargin;     % figure size on paper (widht & hieght)

set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[x y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')

print -dpdf -r0 pHistData.pdf


%=========================================================================
%           FIGURE 2: LARGEST EIGENVALUE (Companion Form)
%=========================================================================

pnames = strvcat('Largest Eigenvalue (Recursive Average)',...
    'Largest Eigenvalue (Posterior Marginal Density)');

figure('Position',[20,20,900,600],'Name',...
    'Largest Eigenvalue (Companion Form)','Color','w')

subplot(1,2,1) 
plot(rmean,'LineStyle','-','Color','k','LineWidth',2.5)
hold on
plot(r5per,'LineStyle','--','Color','k','LineWidth', 0.75)
hold on
plot(r95per,'LineStyle','--','Color','k','LineWidth', 0.75)
grid on
title(pnames(1,:),'FontSize',10,'FontWeight','bold');


subplot(1,2,2)
plot(x_den,density,'LineStyle','-','Color','k','LineWidth',2.5)
hold on
xline(lb_eig,'LineStyle','--','Color','k','LineWidth', 0.75)
hold on
xline(ub_eig,'LineStyle','--','Color','k','LineWidth', 0.75)
grid on
title(pnames(2,:),'FontSize',10,'FontWeight','bold');

x = 29.7;                  % A4 paper size
y = 21.0;                  % A4 paper size
xMargin = 1;               % left/right margins from page borders
yMargin = 1;               % bottom/top margins from page borders
xSize = x - 2*xMargin;     % figure size on paper (widht & hieght)
ySize = y - 2*yMargin;     % figure size on paper (widht & hieght)

set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[x y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')

print -dpdf -r0 pEigen.pdf