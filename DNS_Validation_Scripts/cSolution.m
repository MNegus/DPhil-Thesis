function cs = cSolution(ds, Js, epsilon)
%%cSolution
% Solves for the time-dependent quantity c, used in the analytical solution
% for the force contribution from the inner region

    %% Set-up numerical solving
    options = optimoptions('fsolve', 'TolFun', 1e-10, 'TolX', 1e-10);
    cGuess = sqrt(pi * ds ./ Js) / epsilon;
    zero_fun = @(c) c.^2 + 4 * c + 2 * log(c) + 1 - pi * ds ./ (epsilon^2 * Js);
    
    %% Numerically solves for cs
    cs = fsolve(zero_fun, cGuess, options);

end