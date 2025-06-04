function [SigmaStar,LambdaStar,SigmaIterMat,LambdaIterMat,objValConv,objValStar,objValAvgConv,SigmaAvgIterMat,LambdaAvgIterMat,tdVec]= ...
    FM_MaxMin(LagrangianEq,gradLambda,SigmaHat,epsD,LMOcase,numIter) 
% Solving Max-Min problem using GA with 1/sqrt(t) step size         
        if LMOcase == 1
            fprintf('Running FM_MaxMin: \t with Frobenius norm: \n')
        elseif LMOcase == 2
            fprintf('Running FM_MaxMin: \t with Kullback-Leibler divergence: \n')
        elseif LMOcase == 3
            fprintf('Running FM_MaxMin: \t with Gelbrich distance: \n')
        end

        iterI= 1; % The index of iteration (Iteration number)

        % ****** Constants
        n_xi= size(SigmaHat,1);

        % ****** Dual decision variable's (Lambda's) initial value
        Lambda= generateInitialLambda(n_xi);
        Lambda= (Lambda+Lambda')/2; 
        Lambdasave= -eye(n_xi,n_xi);

        % ****** Max-Min Optimization 
        if iterI == 1
            nonSymmetricLambdaTimes= 0;
        end

        while iterI <= numIter
            if rem(iterI,1000) == 0
                fprintf('Iteration: %d \t',iterI)
                if rem(iterI,5000) == 0
                    fprintf('\n')
                end
            end
           
            if norm(Lambda-Lambda') >= 10^-8 % Checking if Lambda is nonsymmetric
                nonSymmetricLambdaTimes= nonSymmetricLambdaTimes + 1; 
            end

            % ****** Updating Sigma  
            if LMOcase == 1 % Frobenius norm
                Sigma= oracleFrobenius(SigmaHat,Lambda,epsD); % This oracle has been written for min. problem including Lambda
            end
            if LMOcase == 2 % Kullback-Leibler divergence
                Sigma= oracleKullbackLeibler(SigmaHat,-Lambda,epsD); % This oracle has been written for max. problem including D
            end
            if LMOcase == 3 % Gelbrich distance
                Sigma= oracleWassDistRobustKF(SigmaHat,-Lambda,epsD);  % This oracle has been written for max. problem including D
            end

            % ****** Alternatively Projected Gradient Ascent 
            % ********* Gradient Ascent            
            deltaLambda= gradLambda(Sigma);
            
            % ************ Step size ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            td= 1/sqrt(iterI);
            tdVec(iterI)= td;

            % ************ Saving the value of Lambda before updating it
            Lambdasave= Lambda;

            % ************ Updating Lambda ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            LambdaPlus= Lambda + td * deltaLambda;
            
            % *************** Alternating Projection
            Lambda= DykstraProjection(LambdaPlus); 

            % ************ Computing LambdaAvg
            iterAvg= iterI;
            if iterAvg == 1
                LambdaAvg= Lambdasave;
            else
                LambdaAvg= (iterAvg-1)/iterAvg*LambdaAvg + 1/iterAvg*Lambdasave;
            end
            % ************ Computing SigmaAvg
            if iterAvg == 1
                SigmaAvg= Sigma;
            else
                SigmaAvg= (iterAvg-1)/iterAvg*SigmaAvg + 1/iterAvg*Sigma;
            end
            % ************ Saving iterative values of Sigma, Lambda, LambdaAvg, and Corresponding Objective functions  
            SigmaIterMat(:,:,iterI)= Sigma;
            LambdaIterMat(:,:,iterI)= Lambdasave;

            SigmaAvgIterMat(:,:,iterI)= SigmaAvg;
            LambdaAvgIterMat(:,:,iterI)= LambdaAvg;

            objValConv(iterI)= LagrangianEq(Sigma,Lambdasave);
            objValAvgConv(iterI)= LagrangianEq(SigmaAvg,LambdaAvg); 

            % ************ Checking stoping condition
            if iterI >= 2
                if abs((objValConv(iterI) - objValConv(iterI-1)) / objValConv(iterI)) <= 10^-6
                    iterI
                    break
                end
            end
            iterI = iterI + 1;
        end
        
        SigmaStar= Sigma;
        LambdaStar= Lambdasave;
        
        objValStar= LagrangianEq(SigmaStar,LambdaStar);
        fprintf('\n')
end
      


