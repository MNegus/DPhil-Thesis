function [xsLower, hsLower, xsUpper, hsUpper] = outerouterfreesurface(t, SubstrateFunctions)
%OUTERFREESURFACE Free surface location in the outer-outer region
%   Detailed explanation goes here

    %% Load function parameters
    d = SubstrateFunctions.d(t); % Turnover point
    epsilon = SubstrateFunctions.epsilon;
    
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional free-surface
        error("Two-dimensional outer-outer not currently supported.");
        
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric free-surface
        
        %% Clusters lower points near the turnover point and the droplet radius
        xMin = epsilon * d;
        xMax = 1;
        lambda = 10;
        sigmas = linspace(-1, 1, 1e3);
        a = 0.5 * (xMin + xMax);
        b = 0.5 * coth(lambda) * (xMax - xMin);

        xsLower = a + b * tanh(lambda * sigmas);
        
        %% Cluster upper points near the droplet radius
        xMin = 0;
        xMax = 1;
        lambda = 10;
        sigmas = linspace(1, 0, 1e3);
        b = (xMax - xMin) / (exp(-lambda) - 1);
        a = xMin - b;

        xsUpper = a + b * exp(-lambda * sigmas);
        
        %% Find outer-outer free surface
        G = 4 * d^5 / 45;
        Hminus = - G ./ (2 * sqrt(2) * sqrt(1 - xsLower.^2) .* (1 - sqrt(1 - xsLower.^2)).^(3/2));
        Hplus = G ./ (2 * sqrt(2) * sqrt(1 - xsUpper.^2) .* (1 + sqrt(1 - xsUpper.^2)).^(3/2));
        
        % Adjust for singularity
        Hminus(end) = 0;
        Hplus(1) = 0;
        
        %% Original variable solution
        hsLower = 1 - sqrt(1 - xsLower.^2) - epsilon^2 * t + epsilon^5 * Hminus;
        hsUpper = 1 + sqrt(1 - xsUpper.^2) - epsilon^2 * t + epsilon^5 * Hplus;
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'"); 
    end
    
end