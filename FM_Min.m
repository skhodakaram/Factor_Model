function [SigmaStar,LStar,DStar,optObjVal]= FM_Min(SigmaHat,epsD) 
% Min. problem with LMI constraint corresponding to Gelbrich distance - solved by MOSEK
    fprintf('Running FM_Min: \n')
    n_xi= size(SigmaHat,1);
    L= sdpvar(n_xi,n_xi);
    S= sdpvar(n_xi,n_xi);
    C= sdpvar(n_xi,n_xi,'full');
    distLMI= [S  C;
              C' SigmaHat];

    objFunc = trace(L*eye(n_xi));
    Constraint = [L >= 0 , S-L >= 0 , S-L - diag(diag(S-L)) == 0, trace((S + SigmaHat - 2*C)'*eye(n_xi)) <= epsD^2 , distLMI >= 0];

    ops= sdpsettings('solver','mosek','verbose',0);
    optimize(Constraint,objFunc,ops);

    LStar= value(L);
    SigmaStar= value(S);
    CStar= value(C);

    DStar= SigmaStar - LStar;
    optObjVal= trace(LStar*eye(n_xi));
end