function ps = outerpressure(rs, d, d_t, d_tt, epsilon)
%%outerpressure
% Returns the outer solution for the pressure along the substrate in the
% axisymmetric, plate impact case.

    % Define outer pressure to be zero past the turnover point
    ps = zeros(size(rs));
    ps(rs >= epsilon * d) = 0;
    
    % Outer x variable
    rhats = rs(rs < epsilon * d) / epsilon;
    
    % Outer pressure
    ps(rhats < d) = (4 / (3 * pi * epsilon)) ...
        * (d_t^2 * (2 * d^2 - rhats.^2) + d * d_tt * (d^2 - rhats.^2)) ...
        ./ sqrt(d^2 - rhats.^2);
end