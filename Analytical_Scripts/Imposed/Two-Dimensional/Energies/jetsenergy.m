function Es = jetsenergy(ts, SubstrateFunctions)
%JETSENERGY Energy in the jets (splash sheet in axisymmetric)
    
    % Load in function values
    epsilon = SubstrateFunctions.epsilon;
    
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional energy
        % Load substrate functions
        Bs = SubstrateFunctions.B(ts);
        Cs = SubstrateFunctions.C(ts);

        % Define the integrand
        integrand = (epsilon^2 * pi / 2) * (1 - Bs) .* Cs;
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric energy
        % Load substrate functions
        ds = SubstrateFunctions.d(ts);
        d_ts = SubstrateFunctions.d_t(ts);
        
        % Define integrand 
        integrand = (8 * epsilon^3 / 9) * ds.^4 .* d_ts.^3;
        
        % At t == 0, set the integrand to be 0
        zeroIdx = find(ts == 0);
        if ~isempty(zeroIdx)
           integrand(zeroIdx) = 0;
        end
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'");
    end
    
    % Numerically integrate
    Es = cumtrapz(ts, integrand);
    
end

