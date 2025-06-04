clear all
close all

tStart = tic; 

%% *** Flags (Code control parameters): problemCase, LMOcase & SweetSpot 
% ****** Determining the optimiztion problem to be solved
problemCase= 2; % Choices: 1, 2, 3, or 4
                % 1: FM_Min: Min. problem with LMI constraint corresponding to Gelbrich distance - solved by MOSEK
                % 2: FM_MaxMin: Solving Max-Min problem using GA (Gradient Ascent) with 1/sqrt(t) step size
                                LMOcase= 3; % Choices: 1, 2, or 3 
                                            % 1: Frobenius norm
                                            % 2: Kullback-Leibler divergence
                                            % 3: Gelbrich distance
                % 3: FM_Min_Fro: Min. with Frobenius norm constraint - solved by MOSEK
                % 4: FM_Min_KL: Min. with Kullback-Leibler divergence constraint - solved by MOSEK
% ****** NO sweet spot diagram or sweet spot diagram
SweetSpot= 0; % Either 0 or 1 
              % 0: NO sweet spot diagram
              % 1: We plot the sweet spot diagram (We do 'NExp' experiments, in each of which we get 
              %    'N' measurement samples and solve the problem for each element of 'epsMat'.

%% *** Dimension of xi 
n_xi= 50 

%% *** Constants 
% ****** Total number of iterations in FM_MaxMin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
numIter= 10^6;
% ****** Epsilon values (the radius of ball)
epsMat= [0.01, 0.01*sqrt(10), 0.1, 0.1*sqrt(10), 1, 1*sqrt(10), 10, 10*sqrt(10), 100, 100*sqrt(10), 1000];
% ****** Number of Experiments & indice of epsMat
if SweetSpot == 0
    NExp= 1;
    jStart= 7;
    jEnd= jStart;
end
if SweetSpot == 1
    NExp= 100; 
    jStart= 1;
    jEnd= size(epsMat,2);
end

%% *** Objective function 
% ****** Objective function definition
LagrangianEq= @(SigmaEq,LambdaEq) trace(LambdaEq' * SigmaEq);
% ****** Gradient definition
gradLambda= @(SigmaEq) SigmaEq';

%% *** MAIN 
for i= 1:NExp 
    disp(['Experiment:', num2str(i)])
    % ****** Generating Samples
    [xiMat,SigmaTrue,LTrue,DTrue,N]= genSamples(SweetSpot,n_xi);
    % ****** Computing 'SigmaHat' matrix
    SigmaHat= computeSigmaHat(n_xi,xiMat,N);
    % *** Distance ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KL_Div= @(X,XHat) 0.5 * ( -log(det(X)) + log(det(XHat)) + trace(X*inv(XHat)) - n_xi );
    Gel_Dis= @(X,XHat) ( trace(X) + trace(XHat) - 2*trace( (XHat ^ 0.5 * X * XHat^0.5)^0.5 ) )^0.5;
    
    for j = jStart:jEnd
        epsD= epsMat(j);
        fprintf('Optimization for eps= %d \t',epsD)
        % ****** Optimization 
        % ********* Min problem with Gelbrich distance constraint
        if problemCase == 1
            [SigmaStar,LStar,DStar,optObjValMat(i,j)]= FM_Min(SigmaHat,epsD);
        end
        % ********* Max-Min problem
        clear objValConvMat objValAvgConvMat
        if problemCase == 2
            if epsD == 0 % For the case when eps=0, the ball contains only SigmaHat
                SigmaStar= SigmaHat; 
                optObjValMat(i,j)= 0; % We do not have Lambda in this case 
            else
                [SigmaStar,LambdaStar,SigmaIterMat,LambdaIterMat,objValConvMat(i,j,:),optObjValMat(i,j),objValAvgConvMat(i,j,:),SigmaAvgIterMat,...
                    LambdaAvgIterMat,tdVec]= FM_MaxMin(LagrangianEq,gradLambda,SigmaHat,epsD,LMOcase,numIter);
            end
        end
        % ********* Min problem with Frobenius norm constraint 
        if problemCase == 3
            [SigmaStar,LStar,DStar,optObjValMat(i,j)]= FM_Min_Fro(SigmaHat,epsD); 
        end
        
        % ********* Min problem with Kullback-Leibler divergence constraint 
        if problemCase == 4
            [SigmaStar,LStar,DStar,optObjValMat(i,j)]= FM_Min_KL(SigmaHat,epsD); 
        end
        
        % ****** Calculation of the estimation error of SigmaStar and SigmaHat
        [SErrorMat(i,j),SHatErrorMat(i,j)]= calErrorMat(LMOcase,KL_Div,Gel_Dis,SigmaTrue,SigmaStar,SigmaHat);
    end
end

%% *** Plot 
% ****** Plotting sweet spot diagrams
if SweetSpot == 1 
    plotSweetSpot(epsMat,SErrorMat,SHatErrorMat,LMOcase)
end

if SweetSpot == 0
    disp('Optimal objective value:')            
    disp(optObjValMat)
end

tEnd= toc(tStart)

save var.mat  

