function [xiMat,SigmaTrue,LTrue,DTrue] = genSampleMultiNormDist(n_xi,r,N,SweetSpot)
% Generates samples using covariance matrix of 'alpha', which is identity 
% (Covariance matrix of the factors), and covariance matrix of 'omega' (noise) 
% and matrix 'Phi'.
    
    % *** Generating samples of alpha ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mu_alpha = zeros(r,1);
    rng(1,'twister');
    Sigma_alpha = eye(r); % Covariance matrix of alpha
    Phi_True= 5 + rand(n_xi,r); 
    v1= rand(r,1);
    if SweetSpot == 1
        rng('shuffle')
    end
    alphaMat = (mvnrnd(mu_alpha, Sigma_alpha, 2*N))'; % Generates N samples for alpha with mean: mu_alpha and covariance matrix: Sigma_alpha
    
    % *** Generating samples of omega ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    mu_omega = zeros(n_xi,1);
    rng(0);
    DTrue= diag(5 + rand(n_xi,1));
    if SweetSpot == 1
        rng('shuffle')
    end
    omegaMat= (mvnrnd(mu_omega, DTrue, 2*N))';
    
    % *** Constructing matrices L and Sigma ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    LTrue= Phi_True * Sigma_alpha * Phi_True'; % The true L matrix
    SigmaTrue= LTrue + DTrue; % The true Sigma matrix

    % *** Constructing synthetic data matrix ~~~~~~~~~~~~~~~~~~~~~~
    xiTOMat= Phi_True*alphaMat + omegaMat;
    xiMat= xiTOMat(:,1:N);
end

    
  