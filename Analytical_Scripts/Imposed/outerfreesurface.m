function hHats = outerfreesurface(xHats, t, SubstrateFunctions)
%OUTERFREESURFACE Summary of this function goes here
%   Detailed explanation goes here

    %% Load function values
    d = SubstrateFunctions.d(t);
    
    if SubstrateFunctions.dimension == "2D"
        %% Two-dimensional free-surface
        
        % Load function values
        a = SubstrateFunctions.a(t);
        b = SubstrateFunctions.b(t);
        
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
        hHats = 0.5 * xHats .* sqrt(xHats.^2 - d^2) ...
            + sqrt(xHats.^2 - d^2) .* (a * I2_0 + (b / 3) * I2_1) ...
            - (xHats ./ sqrt(xHats.^2 - d^2)) ...
                .* (t - d^2 / 4 - a * I1_0 - (b / 3) * I1_1);
        
        %% Setting turnover point value (Wagner condition)
        dIdx = find(xHats == d);
        if ~isempty(dIdx)
           hHats(dIdx) = -a - b * d^2;
           hHats(dIdx)
        end
        
    elseif SubstrateFunctions.dimension == "axi"
        %% Axisymmetric free-surface
        % Load function values
        w = SubstrateFunctions.w(t);
        
        %% Free surface solution
        hHats = 0.5 * xHats.^2 - t + (d / pi) * sqrt(xHats.^2 - d^2) ...
            + (2 / pi) * asin(d ./ xHats) .* (t - w - 0.5 * xHats.^2);
        
        %% Setting turnover point value (Wagner condition)
        dIdx = find(xHats == d);
        if ~isempty(dIdx)
           hHats(dIdx) = -w;
        end
    else
        error("Invalid dimension. Needs to either be '2D' or 'axi'"); 
    end
        
    % Setting all values before the turnover point to nan
    hHats(xHats < d) = nan;
    
end

