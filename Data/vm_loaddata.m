%Data set:
%Output, Inflation, Interest Rates, InvVelocity  

load datampshock.txt
YY=datampshock;
ti=linspace(1964,0.2, size(YY,1))';
nobs=size(YY,1);

%Convert velocity into real money balances

YY(:,4)=YY(:,4)+YY(:,1);


%% Graphs

% Attributes 
linesty={'-','-','-','-'};
color={'k','b','r','g'};
marker={'none','none','none','none'};
legend_list = {'Output', 'Inflation', 'Interest rate', 'Inverse velocity'};

% xaxis
time = 1965:0.25:2006;

figure(1);
p = plot(time, YY, 'linewidth', 1);
size_p=size(p);
grid on;
for i=1:size_p(1);
    set(p(i),'LineStyle',linesty{i},'Color',color{i},'Marker',marker{i});
    xlabel("Time")
end
h=legend(legend_list, 'orientation', 'horizontal');
set(h,'Fontsize',10);
set(h,'Position',[0.45 0.01 0.1 0.05]);

X = 21.0;                  % A4 paper size
Y = 29.7;                  % A4 paper size
xMargin = 1;               % left/right margins from page borders
yMargin = 1;               % bottom/top margins from page borders
xSize = X - 2*xMargin;     % figure size on paper (widht & hieght)
ySize = Y - 2*yMargin;     % figure size on paper (widht & hieght)

set(gcf, 'Units','centimeters', 'Position',[0 0 xSize ySize]/2)

set(gcf, 'PaperUnits','centimeters')
set(gcf, 'PaperSize',[X Y])
set(gcf, 'PaperPosition',[xMargin yMargin xSize ySize])
set(gcf, 'PaperOrientation','portrait')

print -dpdf -r0 pRawData.pdf