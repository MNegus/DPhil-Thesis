function Es = outerenergy(ts, SubstrateFunctions)
%OUTERENERGY Energy in the outer region
    
    % Load in substrate coefficients
    epsilon = SubstrateFunctions.epsilon;
    ds = SubstrateFunctions.d(ts);
    
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional energy
        % Load substrate coefficients
        a_ts = SubstrateFunctions.a_t(ts);
        b_ts = SubstrateFunctions.b_t(ts);
        
        % Return outer energy
        Es = (pi * epsilon^2 * ds.^2 / 4) .* (1 - a_ts).^2 ...
            + (5 * pi * epsilon^2 * ds.^4 / 192) .* b_ts .* (6 * a_ts + ds.^2 .* b_ts - 6);
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric energy
        % Load substrate coefficients
        w_ts = SubstrateFunctions.w_t(ts);
        d_ts = SubstrateFunctions.d_t(ts);
        
        % Determine energy
        Es = (4 * epsilon^3 / 9) * ds.^4 .* d_ts .* (1 - w_ts);
        Es(1) = 0;
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'");
    end
end

