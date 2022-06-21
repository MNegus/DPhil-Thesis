function hs = outerfreesurface(rs, t, d, SubstrateCoefficients, epsilon)
%OUTERFREESURFACE Summary of this function goes here
%   Detailed explanation goes here

    %% Load substrate coefficients
    w = SubstrateCoefficients.ws;
    
    %% Define xHats
    rHats = rs / epsilon;
    
    %% Free surface solution
    hs = 0.5 * rHats.^2 - t + (d / pi) * sqrt(rHats.^2 - d^2) ...
        + (2 / pi) * asin(d ./ rHats) .* (t - w - 0.5 * rHats.^2);
    
    %% Cutting off at turnover point
    
    % Setting turnover point value (Wagner condition)
    dIdx = find(rHats == d);
    if ~isempty(dIdx)
       hs(dIdx) = -w;
    end
    
    % Setting all values before the turnover point to nan
    hs(rHats < d) = nan;
    
    % Scaling
    hs = epsilon^2 * hs; 
end

