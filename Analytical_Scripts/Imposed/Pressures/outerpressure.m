function ps = outerpressure(xs, d, A, aHat_tt, bHat_tt, epsilon)
%%outerpressure
% Returns the outer solution for the pressure along the substrate in the
% two-dimensional, quadratic substrate case.

    % Define outer pressure to be zero past the turnover point
    ps = zeros(size(xs));
    ps(xs >= epsilon * d) = 0;
    
    % Outer x variable
    xhats = xs(xs < epsilon * d) / epsilon;
    
    % Plate term
    plateTerm = 0.5 * aHat_tt * (d^2 - 2 * xhats.^2);
    
    % Quadratic term
    quadTerm = (1 / 24) * bHat_tt * (d^4 + 4 * d^2 * xhats.^2 - 8 * xhats.^4);
    
    ps(xhats < d) = (1 / epsilon) * (A - plateTerm - quadTerm) ...
        ./ sqrt(d^2 - xhats.^2);
end