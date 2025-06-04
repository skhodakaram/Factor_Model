function Lambda = projection_to_S1(Lambda)
    d = diag(Lambda);
    d(d <= 0) = 0;
    Lambda = Lambda - diag(d);
end