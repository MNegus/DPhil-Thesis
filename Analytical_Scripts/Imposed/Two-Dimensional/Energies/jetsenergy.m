function Es = jetsenergy(ts, SubstrateFunctions, epsilon)
%OUTERENERGY Energy in the outer region
    
    % Load in function values
    Bs = SubstrateFunctions.B(ts);
    Cs = SubstrateFunctions.C(ts);

    % Define the integrand function
    integrand = (1 - Bs) .* Cs;
    
    % Numerically integrate
    Es = (epsilon^2 * pi / 2) * cumtrapz(ts, integrand);
    
end

