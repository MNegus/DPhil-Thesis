function ps = overlappressure(xs, t, SubstrateFunctions)
%%outerpressure
% Returns the overlap solution for the pressure along the substrate in the
% two-dimensional, plate impact case. 

    % Load in parameters
    d = SubstrateFunctions.d(t);
    C = SubstrateFunctions.C(t);
    epsilon = SubstrateFunctions.epsilon;

    % Define overlap pressure to be zero past the turnover point
    ps = zeros(size(xs));
    ps(xs >= epsilon * d) = 0;
    
    % Outer x variable
    xHats = xs(xs < epsilon * d) / epsilon;
    
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional pressure
        ps(xHats < d) = C ./ (epsilon * sqrt(2 * d * (d - xHats)));
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric pressure
        ps(xHats < d) = C ./ (epsilon * sqrt(d - xHats));      
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'");
    end
end