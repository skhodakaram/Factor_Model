function [eigVec,eigVal]= eVec_eVal_sort_Descending(Matrix)
% The eigenvalues of 'Matrix' are sorted in decreasing order in the 'eigVal'. 
% The corresponding eigenvectors are also sorted in 'eigVec'
    % *** Eigendecomposition of Matrix  
    [eigVec,eigVal]= eig(Matrix); 
    % *** Checking if the matrix is nonsymmetric
    if norm(imag(eigVal)) > 10^-8 % We stop running this code if the matrix is not symmetric
        fprintf('Error: The matrix is not (real) symmetric. \n')
        return
    end
    % *** Sorting eigenvalues in decreasing order and their corresponding eigenvectors 
    n= size(Matrix,1);
    for i= 1:n
        max= eigVal(i,i);
        for j= i+1:n
            if eigVal(j,j) > max
                max= eigVal(j,j);
                eValSave= eigVal(i,i);
                eigVal(i,i)= eigVal(j,j);
                eigVal(j,j)= eValSave;
                eVecSave= eigVec(:,i);
                eigVec(:,i)= eigVec(:,j);
                eigVec(:,j)= eVecSave;
            end
        end
    end
end