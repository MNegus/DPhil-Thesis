function ps = overlappressure(xs, d, C, epsilon)
%%outerpressure
% Returns the overlap solution for the pressure along the substrate in the
% two-dimensional, plate impact case. 

    % Define overlap pressure to be zero past the turnover point
    ps = zeros(size(xs));
    ps(xs >= epsilon * d) = 0;
    
    % Overlap pressure on the contact set
    xhats = xs(xs < epsilon * d) / epsilon;
    ps(xhats < d) = C ./ (epsilon * sqrt(2 * d * (d - xhats)));
end