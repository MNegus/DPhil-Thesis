function p = outerpressureglobal(rs, z, d, d_t, d_tt)
%OUTERPRESSUREGLOBAL
%   Assume a constant value of z and an array of rs

    %% Function for the coefficient function, nu_tt
    nu_tt = @(lambda) (4 ./ (3 * pi * lambda.^2)) ...
        .* (d * (d_t^2 + d * d_tt) * lambda .* cos(lambda * d) ...
            - sin(lambda * d) .* ((1 + d^2 * lambda.^2) * d_t^2 + d * d_tt));
        
    %% Function for the integrand
    integrand = @(lambda) ...
        nu_tt(lambda) .* exp(-z * lambda) .* besselj(0, lambda * rs);

    %% Determine p using numerical integration
    lambdas = exp(linspace(-20, 20, 1e3))';
    p = -trapz(lambdas, integrand(lambdas), 1);
%     p = integral(integrand, 0, Inf);
end