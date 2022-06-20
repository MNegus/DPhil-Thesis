function Es = jetsenergy(ts, Bs, Cs, epsilon)
%OUTERENERGY Energy in the outer region
    
    % Define the integrand function
    integrand = (1 - Bs) .* Cs;
    
    % Numerically integrate
    Es = (epsilon^2 * pi / 2) * cumtrapz(ts, integrand);
    
end

