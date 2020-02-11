%==========================================================================
%           ECONOMETRICS IV - SPRING 2020 - Prof. Schorfehide 
%                       REDUCED-FORM VARs                     
%
%
% Author: Juan Castellanos Silvan (from Luigi Boccola's original code)
% Date  : 02/11/2020
%==========================================================================

% Exercise 6: point and interval forecasts

%=========================================================================
%                             HOUSEKEEPING
%=========================================================================

tic
close all
clear all
clc

%=========================================================================
%       DIRECT MONTE CARLO SAMPLING FROM POSTERIOR OF VAR PARAMETERS
%=========================================================================

load("resultsVAR.mat", 'n','p','YYact','Phip','Sigmap')


%=========================================================================
%                      FORECASTING ALGORITHM
%=========================================================================
% Pre-allocation
h = 12;                 % forcasting horizon
nsim = 1000;            % number of simulations
YYf = zeros(h, n, nsim);  
XXf = zeros(h, n*p+1, nsim);

for isim = 1:nsim
    
    % -------------------
    % Initial conditions
    % -------------------
    j=1;
    for i = 1:p
        XXf(1,j:i*p, isim) = YYact(end-i+1,:);
        j = i*p+1;
    end
    XXf(1,end, isim) = 1;
    
    % --------------------
    % Multi-step forecast
    % --------------------
    U = zeros(h,n);
    for ih = 1:h
        
        % Updating
        if ih > 1
            XXf(ih,1:p, isim) = YYf(ih-1,:,isim);
            XXf(ih,p+1:end-1,isim) = XXf(ih-1,1:end-p-1, isim);
            XXf(ih,end, isim) = 1;
        end
        
        % Independent draws for error term
        u = mvnrnd(zeros(n,1), squeeze(Sigmap(isim,:,:)));
        U(ih,:) = u; 
        
        % Reduced-From VAR Model
        YYf(:,:,isim) = XXf(:,:,isim)*squeeze(Phip(isim,:,:)) + U;
    
    end
    
    
end

YY = [repmat(YYact,1,1,nsim); YYf];

%% Fan Charts

figure_list = {'Output', 'Inflation', 'InterestRate', 'InverseVelocity'};

% xaxis 
time = 1966:0.25:2009;

for i=1:n
    figure('Name',figure_list{i})
    [lh, ph] = fanChart(time, squeeze(YY(:,i,:)), 'median');
    set(lh,'LineWidth',1);
    xlabel("Time")
    grid on
    
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
    
   saveas(gcf, strcat('p',figure_list{i},'.pdf'))

end

disp(['         ELAPSED TIME:   ', num2str(toc)]);

elapsedtime=toc;

