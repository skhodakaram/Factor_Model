function []= plotErrorObjVal(objValConvMat,i,j)
    maxIter= size(objValConvMat,3); % Total number of iterations
    tConv= 1:maxIter;

    % *** Plottig objective value
    for c= 1:maxIter
        objValMatPlot(c)= objValConvMat(i,j,c);
    end

    optObjVal= objValMatPlot(maxIter);

    convObjValMatPlot= abs(objValMatPlot-optObjVal*ones(1,maxIter))./abs(optObjVal);
    figure()
        loglog(tConv, convObjValMatPlot)
        set(gca,'FontSize',12)
        xlabel('Iteration($t$)','interpreter','latex','FontSize',14)
        ylabel('$\frac{| \langle \Lambda_t,\Sigma_t \rangle - \langle \Lambda^\star,\Sigma^\star \rangle |}{| \langle \Lambda^\star,\Sigma^\star \rangle |}$','interpreter','latex','FontSize',20)
        grid on
end
