function hs = innerfreesurface_upper(xs, t, SubstrateFunctions)
%%INNERFREESURFACE_UPPER Upper side of the free surface in inner region

    % Load in parameters
    d = SubstrateFunctions.d(t);
    d_t = SubstrateFunctions.d_t(t);
    J = SubstrateFunctions.J(t);
    epsilon = SubstrateFunctions.epsilon;
    
    if SubstrateFunctions.dimension == "2D"
        w = SubstrateFunctions.a(t) + SubstrateFunctions.b(t) * d^2;
    else
        w = SubstrateFunctions.w(t);
    end

    %% Define inner variables
    turnoverIdx = sum(xs < epsilon * d); % Index of the turnover point
    xTildes = (xs - epsilon * d) / epsilon^3;
    
    %% Determine xi paramters
    % The xi parameters for the lower free surface start at xi = -infinity
    % and going up to 0. We need to find the minimum xi value to match
    % with the maximum value of xTilde
    
    xTildeMax = max(xTildes);
    xiMinFun = @(xi) xTildeMax - (J / pi) * (exp(xi) - xi - 1);
    xiMin = fsolve(xiMinFun, - pi * xTildeMax / J);
    
    xis = linspace(xiMin, 0, 3 * length(xTildes));
    
    %% Inner variable free-surface
    xTildesParams = (J / pi) * (exp(xis) - xis - 1);
    hTildes = (J / pi) * (pi + 4 * exp(xis / 2));
    
    %% Original variable free-surface
    hs = -epsilon^2 * w + epsilon^3 * interp1(xTildesParams, hTildes, xTildes);
    
end