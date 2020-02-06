%==========================================================================
%           ECONOMETRICS IV - SPRING 2020 - Prof. Schorfehide 
%                       REDUCED-FORM VARs                     
%
%
% Author: Juan Castellanos Silvan (from Luigi Boccola's original code)
% Date  : 02/04/2020
%==========================================================================

% Question 3: Optimal choice of the hyperparameter

%=========================================================================
%                             HOUSEKEEPING
%=========================================================================

tic
close all
clear all
clc

vmprior = [1 2 3 4 5]';            % indicates the choice of the prior
nAlternatives = length(vmprior);    
vmdd = zeros(nAlternatives, 1);    % stores marginal likelihoods from different tighness

%=========================================================================
%                       LOOP OVER ALTERNATIVES
%=========================================================================
for iAlternative = 1:nAlternatives
   
    %=====================================================================
    %         GENERATE DUMMY OBSERVATIONS FROM MINNESOTA PRIOR 
    %=====================================================================
    disp('                                                                  ');
    disp('      GENERATING DUMMY OBSERVATIONS FROM THE MINESOTA PRIOR       ');
    disp('                                                                  ');
    
    mprior = vmprior(iAlternative);
    vm_dummy

    disp(['      OVERALL TIGHTNESS:   ', num2str(tau)])
    %=====================================================================
    %     DEFINITION OF DATA, LAG STRUCTURE AND POSTERIOR SIMULATION
    %=====================================================================

    [Tdummy,n] = size(YYdum);
    [Tobs,n]   = size(YYact);
    X          = [XXact; XXdum];
    Y          = [YYact; YYdum];
    n          = n;                 % Number of variables in the VAR
    p          = 4;                 % Number of lags in the VAR
    T          = Tobs+Tdummy;
    nsim       = 10000;             % Number of draws from Posterior Density
    nburn      = 0.2*nsim;          % Number of draws to discart
    F          = zeros(n*p,n*p);    % Matrix for Companion Form
    I          = eye(n);

    for i=1:p-1
        F(i*n+1:(i+1)*n,(i-1)*n+1:i*n) = I;
    end


    %=====================================================================
    %               OLS ESTIMATOR FOR PHI AND SSR (SIGMA)
    %=====================================================================

    Phi_tilde = inv(X'*X)*X'*Y;
    Sigma     = (Y-X*Phi_tilde)'*(Y-X*Phi_tilde);

    % Matrices for collecting draws from Posterior Density

    Sigmap    = zeros(nsim,n,n);
    Phip      = zeros(nsim,n*p+1,n);
    largeeig  = zeros(nsim,1);
    counter   = 0;

    %=====================================================================
    %            DRAWS FROM POSTERIOR DENSITY (DIRECT SAMPLING)
    %=====================================================================
    disp('                                                                  ');
    disp('        BAYESIAN ESTIMATION OF VAR: DIRECT SAMPLING...            ');
    disp('                                                                  ');

    for j=1:nsim


        % Draws from the density Sigma | Y

        sigma   = iwishrnd(Sigma,T-n*p-1);

        % Draws from the density vec(Phi) |Sigma(j), Y

        phi_new = mvnrnd(reshape(Phi_tilde,n*(n*p+1),1),kron(sigma,inv(X'*X)));

        % Rearrange vec(Phi) into Phi

        Phi     = reshape(phi_new,n*p+1,n);

        Sigmap(j,:,:) = sigma;
        Phip(j,:,:)   = Phi;

        Phi = Phi(1:n*p,:);

        
         if counter==2000
            disp(['         DRAW NUMBER:   ', num2str(j)]);
            disp('                                                      ');
            disp(['     REMAINING DRAWS:   ', num2str(nsim-j)]);
            disp('                                                      ');

            counter = 0;
         end

    end

    %=====================================================================
    %                        MARGINAL DATA DENSITY
    %=====================================================================

    vm_mdd

    vmdd(iAlternative) = lnpYY;               % Marginal Data Density



    disp(['         ELAPSED TIME:   ', num2str(toc)]);
    disp(' --------------------------------------------')
    

    elapsedtime=toc;

end

% Find the hyperparameter that maximizes likelihood given equal prior prob. 
[mdd_max, iAlternative_max] = max(vmdd);
