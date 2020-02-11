%==========================================================================
%           ECONOMETRICS IV - SPRING 2020 - Prof. Schorfehide 
%                       REDUCED-FORM VARs                     
%
%
% Author: Juan Castellanos Silvan (from Luigi Boccola's original code)
% Date  : 02/06/2020
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

load("resultsVAR.mat")


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

figure_list = {'Output', 'Inflation', 'Interest rate', 'Inverse velocity'};

% xaxis 
time = 1966:0.25:2009;

for i=1:n
    figure('Name',figure_list{i})
    [lh, ph] = fanChart(time, squeeze(YY(:,i,:)), 'median');
    set(lh,'LineWidth',1);
    xlabel("Time")
    grid on
end

disp(['         ELAPSED TIME:   ', num2str(toc)]);

elapsedtime=toc;

