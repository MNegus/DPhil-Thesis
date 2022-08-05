function [xs, hs] = outer_inner_jet_freesurface_composite(xMaxUpper, xMaxLower, t, SubstrateFunctions)
%OUTER_INNER_JET_COMPOSITE Composite solution across outer, inner and jet
%regions of the free surface.
%   Constructs the composite solution of the free surface across the outer,
%   inner and jet regions. Consists of a piecewise solution, where the
%   upper part is the composite between the outer and the upper side of the
%   inner solution, and the lower part is the composite between the inner
%   and the jet region. User needs to specify xMax, which is the maximum
%   value of x chosen. The function will automatically discretise with
%   approximately 1000 points per side.
%   CURRENTLY ONLY SUPPORTED FOR AXISYMMETRIC SOLUTION. THIS COULD BE
%   CHANGED IF THE OUTER-INNER COMPOSITE IS SET TO BE FROM THE INNER
%   SOLUTION, WHICH WILL BE DONE LATER.

    %% Load substrate dependent values
    d = SubstrateFunctions.d(t);
    J = SubstrateFunctions.J(t);
    epsilon = SubstrateFunctions.epsilon;
    w = SubstrateFunctions.w(t);
    
    %% Outer, inner, outer-outer composite
    % Cluster points near turnover point
    xMin = epsilon * d;
    lambda = 10;
    sigmas = linspace(0, 1, 1e3);
    b = (xMaxUpper - xMin) / (exp(lambda) - 1);
    a = xMin - b;

    xsUpper = a + b * exp(lambda * sigmas);
    
    % Outer-Inner composite
    hsOuter = outerfreesurface(xsUpper, t, SubstrateFunctions);
    hsInnerUpper = innerfreesurface_upper(xsUpper, t, SubstrateFunctions);
    hsOverlapOuter = inner_outer_overlapfreesurface(xsUpper, t, SubstrateFunctions);
    
    hsCompositeOuterInner = hsOuter + hsInnerUpper - hsOverlapOuter;
    
    % Outer-outer composite
    hsOuterOuter = outerouterfreesurface(xsUpper, t, SubstrateFunctions);
    hsOverlapOuterOuter = outer_outer_overlap(xsUpper, t, SubstrateFunctions);
    
    hsUpper = hsCompositeOuterInner + hsOuterOuter - hsOverlapOuterOuter;
    
    %% Inner-jet composite
    lambda = 10;
    b = (xMaxLower - xMin) / (exp(lambda) - 1);
    a = xMin - b;
    xsLower = a + b * exp(lambda * sigmas);
    
    hsJet = jetfreesurface(xsLower, t, SubstrateFunctions);
    hsInnerLower = innerfreesurface_lower(xsLower, t, SubstrateFunctions);
    
    % Overlap just set to the jet thickness, needs checking
    hsCompositeJet = hsInnerLower + hsJet - epsilon^3 * J;
    
    %% Add turnover point location as an extra point
    xTurnover = epsilon * d;
    hTurnover = -epsilon^2 * w + (epsilon^3 * J / pi) * (pi + 4);
    
    %% Combine solutions
    % Reverse order of lower solutions
    xs = [flip(xsLower), xTurnover, xsUpper];
    hs = [flip(hsCompositeJet), hTurnover, hsUpper];

end