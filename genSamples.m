function [xiMat,SigmaTrue,LTrue,DTrue,N]= genSamples(SweetSpot,n_xi)
% Generating samples of xi, and generating SigmaTrue, LTrue, and DTrue matrices
    % n_xi is the dimension of the samples
    % ~~~ Constants
    N= 15*n_xi; % the number of samples
    r= ceil(n_xi/5); % the number of factors
    % ~~~ Sample generation
    [xiMat,SigmaTrue,LTrue,DTrue] = genSampleMultiNormDist(n_xi,r,N,SweetSpot);
end