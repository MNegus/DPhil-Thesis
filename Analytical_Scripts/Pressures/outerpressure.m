function ps = outerpressure(xs, t, SubstrateFunctions)
%%outerpressure
% Returns the outer solution for the pressure along the substrate in the
% two-dimensional, quadratic substrate case.

    % Load in parameters
    d = SubstrateFunctions.d(t);
    A = SubstrateFunctions.A(t);
    epsilon = SubstrateFunctions.epsilon;
    
    % Define outer pressure to be zero past the turnover point
    ps = zeros(size(xs));
    ps(xs >= epsilon * d) = 0;
    
    % Outer x variable
    xHats = xs(xs < epsilon * d) / epsilon;
    
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional pressure
        % Load substrate quantities
        a_tt = SubstrateFunctions.a_tt(t);
        b_tt = SubstrateFunctions.b_tt(t);
        
        % Plate term
        plateTerm = 0.5 * a_tt * (d^2 - 2 * xHats.^2);

        % Quadratic term
        quadTerm = (1 / 24) * b_tt * (d^4 + 4 * d^2 * xHats.^2 - 8 * xHats.^4);

        ps(xHats < d) = (1 / epsilon) * (A - plateTerm - quadTerm) ...
            ./ sqrt(d^2 - xHats.^2);
        
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric pressure (xHats -> rHats)
        d_t = SubstrateFunctions.d_t(t);
        d_tt = SubstrateFunctions.d_tt(t);
        
        ps(xHats < d) = (4 / (3 * pi * epsilon)) ...
        * (d_t^2 * (2 * d^2 - xHats.^2) + d * d_tt * (d^2 - xHats.^2)) ...
        ./ sqrt(d^2 - xHats.^2);
    
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'");
    end
end