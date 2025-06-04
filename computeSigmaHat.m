function [SigmaHat]= computeSigmaHat(n_xi,xiMat,N)
% Computing 'SigmaHat' matrix (experimental covariance matrix of 'xi')
%           'SigmaHat' is an estimation for 'SigmaTrue' matrix
    sum= zeros(n_xi,n_xi);
    for i= 1:N
        sum= sum + xiMat(:,i)*xiMat(:,i)';
    end
    SigmaHat= sum / N; 
end