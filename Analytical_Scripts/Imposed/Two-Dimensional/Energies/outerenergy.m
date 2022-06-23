function Es = outerenergy(ts, SubstrateFunctions)
%OUTERENERGY Energy in the outer region
    
    % Load in substrate coefficients
    epsilon = SubstrateFunctions.epsilon;
    ds = SubstrateFunctions.d(ts);
    
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional energy
        % Load substrate coefficients
        aHat_ts = SubstrateFunctions.aHat_t(ts);
        bHat_ts = SubstrateFunctions.bHat_t(ts);
        
        % Return outer energy
        Es = (pi * epsilon^2 * ds.^2 / 4) .* (1 - aHat_ts).^2 ...
            + (5 * pi * epsilon^2 * ds.^4 / 192) .* bHat_ts .* (6 * aHat_ts + ds.^2 .* bHat_ts - 6);
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric energy
        % Load substrate coefficients
        w_ts = SubstrateFunctions.w_t(ts);
        
        % Determine energy
        Es = (4 * epsilon^3 / 9) * ds.^4 .* d_ts .* (1 - w_ts);
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'");
    end
end

