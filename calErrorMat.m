function [SErrorMat,SHatErrorMat]= calErrorMat(LMOcase,KL_Div,Gel_Dis,SigmaTrue,SigmaStar,SigmaHat) 
    if LMOcase == 1 % Frobenius norm
        SErrorMat= norm(SigmaTrue - SigmaStar,'fro');
        SHatErrorMat= norm(SigmaTrue - SigmaHat,'fro');
    end
    if LMOcase == 2 % Kullback-Leibler divergence
        SErrorMat= KL_Div(SigmaStar,SigmaTrue);
        SHatErrorMat= KL_Div(SigmaHat,SigmaTrue);
    end
    if LMOcase == 3 % Gelbrich distance
        SErrorMat= Gel_Dis(SigmaStar,SigmaTrue);
        SHatErrorMat= Gel_Dis(SigmaHat,SigmaTrue);
    end
end