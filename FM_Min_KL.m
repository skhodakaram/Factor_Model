function [SigmaStar,LStar,DStar,optObjVal]= FM_Min_KL(SigmaHat,epsD) 
% Min. with Kullback-Leibler divergence constraint - solved by MOSEK
    fprintf('Running FM_Min_KL: \n')
    n_xi= size(SigmaHat,1);
    L= sdpvar(n_xi,n_xi);
    S= sdpvar(n_xi,n_xi);
    Z= sdpvar(n_xi,n_xi,'full');
    
    for i= 1:n_xi
        for j= i+1:n_xi
            Z(i,j)= 0;
        end
    end
    
    distLMI= [S    Z;
              Z'   diag(diag(Z))];

    sumLogDiagZ= 0;
    for i= 1:n_xi
        sumLogDiagZ= sumLogDiagZ + log(Z(i,i));
    end

    objFunc = trace(L*eye(n_xi)); 
    Constraint = [L >= 0 , S-L >= 0, S-L - diag(diag(S-L)) == 0, distLMI >= 0, log(det(SigmaHat)) + trace(S*inv(SigmaHat)) - n_xi - 2*epsD <= sumLogDiagZ ]; 

    ops= sdpsettings('solver','mosek','verbose',0); 
    optimize(Constraint,objFunc,ops);
    
    LStar= value(L);
    SigmaStar= value(S);
    
    DStar= SigmaStar - LStar;
    optObjVal= trace(LStar*eye(n_xi));
end