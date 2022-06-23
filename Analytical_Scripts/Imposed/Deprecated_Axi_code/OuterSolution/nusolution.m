function [nu, nu_t, nu_tt] = nusolution(lambda, d, d_t, d_tt)
%NUSOLUTION Solution for the coefficient function nu
%   Returns the solution for nu, and its derivatives, which is the
%   coefficient function for the outer region solution in axisymmetric
%   setting. We assume lambda is a vector but d, d_t and d_tt are single
%   values.


    % nu
    nu = (4 ./ (3 * pi * lambda.^4)) ...
        .* ((lambda.^2 * d - 3) .* sin(lambda * d) + 3 * d * lambda .* cos(lambda * d));
    
    % nu_t
    nu_t  = (4 * d * d_t ./ (3 * pi * lambda.^2)) ...
        .* (d * lambda .* cos(lambda * d) - sin(lambda * d));
    
    % nu_tt
    nu_tt = (4 ./ (3 * pi * lambda.^4)) ...
        .* (d * (d_t^2 + d * d_tt) * lambda .* cos(lambda * d) ...
            - sin(lambda * d) .* ((1 + d^2 * lambda.^2) * d_t^2 + d * d_tt));

end