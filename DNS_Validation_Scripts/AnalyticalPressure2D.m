function p0s = AnalyticalPressure2D(ts, ws, w_ts, w_tts, lambda, epsilon)
%%AnalyticalPressure2D
% Determines the analytical pressure at the origin, p0s, in the 2D flat 
% deformable substrate case.

    options = optimoptions('fsolve', 'TolFun', 1e-10, 'TolX', 1e-10);
    
    %% Solution for ds
    dsGuess = 2 * sqrt(ts);
    ds_zero_fun = @(d) ts - d.^2 / 4 - ws .* besselj(0, epsilon * lambda * d);
    ds = fsolve(ds_zero_fun, dsGuess, options);
    
    %% Solgit _ts .* besselj(0, epsilon * lambda * ds)) ... 
    d_ts = (1 - w_ts .* besselj(0, epsilon * lambda * ds)) ...
       ./ (ds / 2 - epsilon * lambda * ws .* besselj(1, epsilon * lambda * ds));
   
    %% Solution for As
    As = ds .* (d_ts .* (1 - w_ts .* besselj(0, epsilon * lambda * ds)) ...
        - w_tts .* besselj(1, epsilon * lambda * ds) / (epsilon * lambda));
    
    %% Determine principal value integral
    integrands = @(xi) sqrt(1 - xi^2) * cos(epsilon * lambda * ds * xi) / xi;
    IntegralSolutions = integral(integrands, -1, 1,'ArrayValued',true);
    
    %% Return pressure
    p0s = (1 ./ ds) .* (As - (ds .* w_tts / pi) .* IntegralSolutions);
    

end