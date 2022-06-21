function Es = outerenergy(ds, d_ts, SubstrateCoefficients, epsilon)
%OUTERENERGY Energy in the outer region
    
    % Load in substrate coefficients
    w_ts = SubstrateCoefficients.w_ts;
    
    % Determine energy
    Es = (4 * epsilon^3 / 9) * ds.^4 .* d_ts .* (1 - w_ts);
end

