function hs = outer_outer_overlap(xs, t, SubstrateFunctions)
%OUTERFREESURFACE Summary of this function goes here
%   Detailed explanation goes here

    %% Load function parameters
    d = SubstrateFunctions.d(t); % Turnover point
    epsilon = SubstrateFunctions.epsilon;
    hs = zeros(size(xs));
    
    
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional free-surface
        error("Two-dimensional outer-outer not currently supported.");
        
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric free-surface
        % Load function values
        G = 4 * d^5 / 45;
        
        %% Original variable overlap
        hs = 0.5 * xs.^2 - epsilon^2 * t - epsilon^5 * G ./ xs.^3;
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'"); 
    end
    
end