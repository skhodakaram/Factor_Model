function [Sigma]= oracleKullbackLeibler(SigmaHat,D,eps) 
% Computing the oracle O(Lambda), when d is defined based on Kullback-Leibler divergence
    
    n_xi= size(SigmaHat,2);
    LB= max(0 , 2 * max( eig(SigmaHat^0.5 * D * SigmaHat^0.5) ));
    UB_1 = norm(svd(- SigmaHat^0.5 * D * SigmaHat^0.5),1);
    if eps <= 1/24
        UB_2= (6/eps)^0.5;
    elseif eps > 1/24
        UB_2= 6 + 1/(4 * eps);
    end
    UB= UB_1 * UB_2; 
    kBisection= 1;
    
    SigmaEq= @(x) inv( inv(SigmaHat) - (2/x)*D );
    KL= @(x) 0.5 * ( -log(det(SigmaEq(x))) + log(det(SigmaHat)) + trace(SigmaEq(x)*inv(SigmaHat)) - n_xi );
    partialObjFunc_gamma= @(x) - KL(x) + eps; 
    
    % *** Bisection method
    while true        
        gamma= (LB+UB)/2;
        
        if UB - LB <= 10^-3 
            break
        end
        if partialObjFunc_gamma(gamma) <= 0
            LB= gamma;
        else
            UB= gamma;
        end
        kBisection= kBisection + 1;
    end

    Sigma= SigmaEq(gamma);
end

 