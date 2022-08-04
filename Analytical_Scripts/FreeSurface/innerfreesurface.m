function [xs, hs, xis] = innerfreesurface(xMax, t, SubstrateFunctions)
%%INNERFREESURFACE Free surface in the inner region

    % Load in parameters
    d = SubstrateFunctions.d(t);
    J = SubstrateFunctions.J(t);
    epsilon = SubstrateFunctions.epsilon;
    
    if SubstrateFunctions.dimension == "2D"
        w = SubstrateFunctions.a(t) + SubstrateFunctions.b(t) * d^2;
    else
        w = SubstrateFunctions.w(t);
    end

    % Define xTilde function
    xTildeMax = (xMax - epsilon * d) / epsilon^3;

    %% Determine xi parameters
%     xiMin = - pi * xTildeMax / J
%     xiMax = pi * xTildeMax / J
    xiMin = -1e3;
    xiMax = 1e3;
    
    xiLowers = linspace(xiMin, 0, 5e2);
    
    xiUppers = linspace(0, xiMax, 5e2);
    xiUppers = xiUppers(2 : end);
    
    xis = [xiLowers, xiUppers];
    
    %% Determine inner free surface location
    xTildesLowers = (J / pi) * (exp(xiLowers) - xiLowers - 1); 
    hTildesLowers = (J / pi) * (pi + 4 * exp(xiLowers / 2));
    
    xTildesUppers = (J / pi) * (xiUppers - log(1 + xiUppers));
    hTildesUppers = (J / pi) * (pi + 4 * sqrt(1 + xiUppers));
    
    
    xTildes = [xTildesLowers, xTildesUppers];
    hTildes = [hTildesLowers, hTildesUppers];
    
    %% Determine solution in original variables
    xs = epsilon * d + epsilon^3 * xTildes;
    hs = -epsilon^2 * w + epsilon^3 * hTildes;
    
    %% Determine sigma parameters
    % The sigma parameters for the free surface go from 0 to infinity, with
    % sigma -> 0 into the jet and sigma -> infinity into the outer. We need
    % to find sigmaMin and sigmaMax to both match with xMax.
    
%     % Solver for finding sigma
%     sigmaZeroFun = @(sigma) xTildeMax - - (J / pi) * (sigma - log(sigma) - 1);
%     
%     % Solve for lower sigma
%     sigmaMinGuess = exp(-pi * xTildeMax / J);
%     sigmaMin = fsolve(sigmaZeroFun, sigmaMinGuess);
%     
%     % Solver for upper sigma
%     
   
    
    %% Inner variable free-surface
%     xTildesParams = (J / pi) * (sigmas - log(sigmas) - 1);
%     hTildes = (J / pi) * (pi + 4 * sqrt(sigmas));
    
    %% Original variable free-surface
%     hs = -epsilon^2 * w + epsilon^3 * interp1(xTildesParams, hTildes, xTildes);
    
end