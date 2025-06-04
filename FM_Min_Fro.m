function [SigmaStar,LStar,DStar,optObjVal]= FM_Min_Fro(SigmaHat,epsD)
% Min. with Frobenius norm constraint - solved by MOSEK
    fprintf('Running FM_Min_Fro: \n')
    n_xi= size(SigmaHat,1);
    L= sdpvar(n_xi,n_xi);
    S= sdpvar(n_xi,n_xi);
    
    objFunc = trace(L*eye(n_xi));
    Constraint = [L >= 0 , S-L >= 0, S-L - diag(diag(S-L)) == 0, norm(S - SigmaHat,'fro') <= epsD];
    
    ops= sdpsettings('solver','mosek','verbose',0); 
    optimize(Constraint,objFunc,ops);
        
    LStar= value(L);
    SigmaStar= value(S);
    
    DStar= SigmaStar - LStar;
    optObjVal= trace(LStar*eye(n_xi));
end