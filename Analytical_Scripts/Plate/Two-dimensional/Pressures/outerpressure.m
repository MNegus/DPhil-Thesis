function ps = outerpressure(xs, d, A, w_tt, epsilon)
%%outerpressure
% Returns the outer solution for the pressure along the substrate in the
% two-dimensional, plate impact case. 

    % Define outer pressure to be zero past the turnover point
    ps = zeros(size(xs));
    ps(xs >= epsilon * d) = 0;
    
    % Outer pressure on the contact set
    xhats = xs(xs < epsilon * d) / epsilon;
    ps(xhats < d) = (1 / epsilon) * (A - 0.5 * w_tt * (d^2 - 2 * xhats.^2)) ...
        ./ sqrt(d^2 - xhats.^2);
end