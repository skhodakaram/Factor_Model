function [Sigma]= oracleFrobenius(SigmaHat,Lambda,epsD)
% Computing the oracle O(Lambda), when d is defined based on Frobenius norm
    % *** Bisection method
    SigmaS_g= @(x) projPosSemiDef(SigmaHat - (1/(2*x))*Lambda);
    derivFunc= @(x) ( norm(SigmaS_g(x) - SigmaHat,'fro') )^2 - epsD^2;
    
    gamma_min= 0;
    gamma_max= norm(Lambda,'fro');
    
    if abs(derivFunc(gamma_max)) <= 10^-10 % Checking if SigmaStar is on the boundary of the Frobenius ball
        gamma= gamma_max;
    else
        % ****** Bisection algorithm for finding gamma*
        LB= gamma_min;
        UB= gamma_max;
        gamma= (LB + UB)/2;
        while UB - LB > 10^-8 
            if derivFunc(gamma) < 0
                UB= gamma;
            else
                LB= gamma;
            end
            gamma= (LB + UB)/2;
        end
    end
    % ****** Coumputing Sigma*(Lambda)
    Sigma= SigmaS_g(gamma);
end

    