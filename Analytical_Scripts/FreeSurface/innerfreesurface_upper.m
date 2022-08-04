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
    
    %% Determine sigma parameters
    % The sigma parameters for the upper free surface start at sigma = 1
    % and go up to infinity. We need to find the maximum sigma value to max
    % with the maximum value of xTilde
    xTildeMax = max(xTildes);
    sigmaMaxFun = @(sigma) xTildeMax - (J / pi) * (sigma - log(sigma) - 1);
    sigmaMax = fsolve(sigmaMaxFun, 1e3);
    
    sigmas = linspace(1, sigmaMax, 3 * length(xTildes));
    
    %% Inner variable free-surface
    xTildesParams = (J / pi) * (sigmas - log(sigmas) - 1);
    hTildes = (J / pi) * (pi + 4 * sqrt(sigmas));
    
    %% Original variable free-surface
    hs = -epsilon^2 * w + epsilon^3 * interp1(xTildesParams, hTildes, xTildes);
    
end