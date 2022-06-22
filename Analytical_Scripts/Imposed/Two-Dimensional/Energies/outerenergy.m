function Es = outerenergy(ts, SubstrateFunctions, epsilon)
%OUTERENERGY Energy in the outer region
    
    % Load in substrate coefficients
    ds = SubstrateFunctions.d(ts);
    aHat_ts = SubstrateFunctions.aHat_t(ts);
    bHat_ts = SubstrateFunctions.bHat_t(ts);
    
    % Return outer energy
    Es = (pi * epsilon^2 * ds.^2 / 4) .* (1 - aHat_ts).^2 ...
        + (5 * pi * epsilon^2 * ds.^4 / 192) .* bHat_ts .* (6 * aHat_ts + ds.^2 .* bHat_ts - 6);
end

