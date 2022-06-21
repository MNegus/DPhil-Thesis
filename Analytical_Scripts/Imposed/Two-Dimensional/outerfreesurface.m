function hs = outerfreesurface(xs, t, d, SubstrateCoefficients, epsilon)
%OUTERFREESURFACE Summary of this function goes here
%   Detailed explanation goes here

    %% Load substrate coefficients
    aHat = SubstrateCoefficients.aHats;
    bHat = SubstrateCoefficients.bHats;
    
    %% Define xHats
    xHats = xs / epsilon;
    
    %% Complex integral solutions
    % I1 integrals
    I1_0 = 1 - xHats ./ (xHats.^2 - d^2).^(1/2);
    I1_1 = d^2 / 2 + xHats.^2 - xHats.^3 ./ (xHats.^2 - d^2).^(1/2);
    
    % I2 integrals
    I2_0 = d^2 ./ (xHats.^2 - d^2).^(3/2);
    I2_1 = xHats .* (2 * xHats.^2 .* (sqrt(xHats.^2 - d^2) - xHats) ...
        + d^2 * (3 * xHats - 2 * sqrt(xHats.^2 - d^2))) ...
        ./ (xHats.^2 - d^2).^(3/2);
    
    %% Free surface solution
    hs = 0.5 * xHats .* sqrt(xHats.^2 - d^2) ...
        + sqrt(xHats.^2 - d^2) .* (aHat * I2_0 + (bHat / 3) * I2_1) ...
        - (xHats ./ sqrt(xHats.^2 - d^2)) ...
            .* (t - d^2 / 4 - aHat * I1_0 - (bHat / 3) * I1_1);
    
    %% Cutting off at turnover point
    
    % Setting turnover point value (Wagner condition)
    dIdx = find(xHats == d);
    if ~isempty(dIdx)
       hs(dIdx) = -aHat - bHat * d^2;
    end
    
    % Setting all values before the turnover point to nan
    hs(xHats < d) = nan;
    
    % Scaling
    hs = epsilon^2 * hs; 
end

