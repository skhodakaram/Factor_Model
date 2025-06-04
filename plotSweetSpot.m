function []= plotSweetSpot(epsMat,SErrorMat,SHatErrorMat,LMOcase)
    % *** Plot tube ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % ****** Constants
    ColorTube1  = [0.0000    0.7    0.7];
    ColorTube4  = [0.3010    0.7450    0.9330];
    Trans = 0.25;
    % ****** Variables
    param.values= epsMat(1:end); 
    param.interval = [5 95];
    % ****** Plot
    figure
        data= SErrorMat(:,:)./SHatErrorMat(:,:); 
        hold on; grid; box on;
        set(gca, 'XScale','log','YScale','log')
        plot(param.values, mean(data,1),'-.','Color',ColorTube4,'LineWidth',2);
        set(gca,'FontSize',11)
        TUBE = PlotTube(data, param, ColorTube1);
        set(TUBE.plot,'facealpha',Trans-0.1);
        legend('$\textrm{mean}$','$[5-95] \textrm{-th percentile}$','interpreter','latex') 
        xlabel('$\varepsilon$','interpreter','latex','FontSize',15)
        if LMOcase == 1 % Frobenius norm
            ylabel('$\frac{\| \Sigma_{True} - \Sigma^* \|_F}{\| \Sigma_{True} - \hat{\Sigma} \|_F}$','interpreter','latex','FontSize',20) 
        end
        if LMOcase == 2 % Kullback-Leibler divergence
            ylabel('$\frac{KL(\Sigma^* || \Sigma_{True})}{KL(\hat{\Sigma} || \Sigma_{True})}$','interpreter','latex','FontSize',20) 
        end
        if LMOcase == 3 % Gelbrich distance
            ylabel('$\frac{G(\Sigma^* , \Sigma_{True})}{G(\hat{\Sigma} , \Sigma_{True})}$','interpreter','latex','FontSize',20) 
        end
end