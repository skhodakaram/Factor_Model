function [MatProj]= projPosSemiDef(Mat) 
% Projection on the positive semidefinite cone 
    [V, D] = eig(Mat);
    d = max(0,diag(D));
    VDh = V*diag(max(0,sqrt(d)));
    MatProj = VDh*VDh';
    return;

    n= size(Mat,1);
    [eigVec,eigVal]= eig(Mat);
    if norm(imag(eigVal)) > 10^-8
        fprintf('Error: The matrix is not (real) symmetric. \n')
        return
    end
    eigVal= eigVal * ones(n,1);
    EigValMatProj= max(eigVal , zeros(n,1));
    if norm(Mat-Mat') < 10^-8
        MatProj= eigVec * diag(EigValMatProj) * eigVec';
    else
        MatProj= eigVec * diag(EigValMatProj) * inv(eigVec);
    end
    if norm(imag(MatProj)) < 10^-10
        MatProj= real(MatProj);
    end
end

