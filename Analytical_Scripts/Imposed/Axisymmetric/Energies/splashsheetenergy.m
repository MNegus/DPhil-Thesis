function Es = splashsheetenergy(ts, ds, d_ts, epsilon)
%OUTERENERGY Energy in the outer region
    
    % Define the integrand function
    integrand = ds.^4 .* d_ts.^3;
    
    % At d == 0, set the integrand to be 0
    zeroIdx = find(ts == 0);
    if ~isempty(zeroIdx)
       integrand(zeroIdx) = 0;
    end
    
    % Numerically integrate
    Es = (8 * epsilon^3 / 9) * cumtrapz(ts, integrand);
    
end

