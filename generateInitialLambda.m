% Author: Shabnam Khodakaramzadeh, researcher at TU Delft
function Lambda= generateInitialLambda(n_xi)
    rng(1,'twister');
    Q= generateSymPosDefMatrix(n_xi);
    Lambda= DykstraProjection(Q);
end