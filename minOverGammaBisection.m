function [gamma]= minOverGammaBisection(partialObjFunc_gamma,LB0,UB0)
% This function is used for finding the minimum of a convex function, objFunc_gamma, using bisection method 
% partialObjFunc_gamma is the partial derivative of objFunc_gamma with respect to gamma

    LB= LB0;
    UB= UB0;
    kBisection= 1;
    while true
        MB= (LB+UB)/2;
        
        if partialObjFunc_gamma(MB) <= 0
            LB= MB;
        else
            UB= MB;
        end
        
        if UB - LB <= 10^-6
            gamma= MB;
            break
        end
        kBisection= kBisection + 1;
    end
end