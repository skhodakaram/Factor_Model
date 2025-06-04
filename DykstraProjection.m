function Lambda = DykstraProjection(Lambda)
    n_xi = size(Lambda, 1);

    Lambda_1 = Lambda;
    Delta_1 = zeros(n_xi, n_xi);
    Delta_2 = zeros(n_xi, n_xi);
    cnt = 0;

    while true
        Lambda_2 = projection_to_S1(Lambda_1 + Delta_1);
        Delta_1 = Lambda_1 + Delta_1 - Lambda_2;
        Lambda_1 = projection_to_S2(Lambda_2 + Delta_2);
        Delta_2 = Lambda_2 + Delta_2 - Lambda_1;
        
        if norm(Lambda_1 - Lambda_2, 'fro') / norm(Lambda_1, 'fro') <= 1e-6
            break
        end
        cnt = cnt + 1;
    end
    
    Lambda = 0.5 * Lambda_1 + 0.5 * Lambda_2;
end



