function Lambda = projection_to_S2(Lambda)
    n= size(Lambda, 1);
    for i= 1 : n
        [v_max, d_max] = eigs(Lambda, 1, 'largestreal', 'Tolerance', 1e-6);
        if isnan(d_max) % d_max == NaN
            [eigVec,eigVal]= eVec_eVal_sort_Descending(Lambda);
            v_max= eigVec(:,1);
            d_max= eigVal(1,1);
        end
        if d_max <= 1 + 1e-6
            break
        else
            Lambda = Lambda + (1 - d_max) * (v_max * v_max');
        end
    end
end