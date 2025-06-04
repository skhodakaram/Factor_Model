function Q= generateSymPosDefMatrix(n)
    eValQ= rand(n,1) + 0.2*ones(n,1); 
    eValQ= diag(eValQ);
    V= rand(n,n);
    for k= 1:n
        sum= zeros(n,1);
        if k == 1
            pV= V(:,k);
        else
            for i= 1:k-1
                sum= sum + (V(:,k)'*eVecQ(:,i))*eVecQ(:,i);
            end
            pV= V(:,k) - sum;
        end
        eVecQ(:,k)= pV / norm(pV,2);
    end

    Q= eVecQ*eValQ*eVecQ';
end