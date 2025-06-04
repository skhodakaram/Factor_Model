function [Sigma]= oracleWassDistRobustKF(SigmaHat,D,epsD) 
% Computing the oracle O(Lambda), when d is defined based on Gelbrich distance

    n_xi= size(SigmaHat,2);
    
    % Eigen decomposition
    [eVec, eVal] = eigs((D+D')/2, length(D), 'la');
    d = diag(eVal)';
    
    L_tilde_gamma= @(x) eVec * diag(1 ./ (x - d)) * eVec';
    partialObjFunc_gamma= @(x) epsD^2 - trace(SigmaHat * (eye(n_xi,n_xi) - x*L_tilde_gamma(x))^2);

    % Bisection
    LB0= max( eig((D+D')/2) ); 
    UB0= norm(-D,'fro') * ( 2 * (max(eig(SigmaHat)))^0.5 + epsD); 

    gamma= minOverGammaBisection(partialObjFunc_gamma,LB0,UB0);

    L_tilde = eVec * diag(1 ./ (gamma - d)) * eVec';

    Sigma= gamma^2 * L_tilde * SigmaHat * L_tilde;
end

    