function ps = overlappressure(rs, d, d_t, epsilon)
%%outerpressure
% Returns the overlap solution for the pressure along the substrate in the
% two-dimensional, plate impact case. 

    % Define overlap pressure to be zero past the turnover point
    ps = zeros(size(rs));
    ps(rs >= epsilon * d) = 0;
    
    % Overlap pressure on the contact set
    rhats = rs(rs < epsilon * d) / epsilon;
    ps(rhats < d) = 2 * sqrt(2) * d^(3/2) * d_t^2 ...
        ./ (3 * pi * epsilon * sqrt(d - rhats));
end