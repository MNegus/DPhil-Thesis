close all;
clear;
addpath("../FreeSurface");

SubstrateFunctions = imposedsubstratefunctions("stationary", "axi");

%% Parameters
ts = linspace(0.001, 1, 1e2);
figure(1);
for t = ts
% for t = 0.5
    
    d = SubstrateFunctions.d(t);
    epsilon = SubstrateFunctions.epsilon;
    J = SubstrateFunctions.J(t);
    xMax = 1.1 * epsilon * d;
    xs = linspace(epsilon * d, xMax, 1e3);

    %% Outer composite
    hsOuter = outerfreesurface(xs, t, SubstrateFunctions);
    hsInnerUpper = innerfreesurface_upper(xs, t, SubstrateFunctions);
    hsOverlapOuter = inner_outer_overlapfreesurface(xs, t, SubstrateFunctions);

    % Composite solution
    hsCompositeOuter = hsOuter + hsInnerUpper - hsOverlapOuter;
    
    %% Jet composite
    [xsJet, hsJet] = jetfreesurface(xMax, t, SubstrateFunctions);
    hsInnerLower = innerfreesurface_lower(xsJet, t, SubstrateFunctions);
    hsOverlapJet = jet_inner_overlapfreesurface(xsJet, t, SubstrateFunctions);
%     hsCompositeJet = hsInnerLower + hsJet - hsOverlapJet;
    hsCompositeJet = hsInnerLower + hsJet - epsilon^3 * J;
    
    %% Plot
    
%     plot(xs, hsOuter, 'linewidth', 2);
%     hold on;
%     plot(xsInner, hsInner, 'linewidth', 2);
%     plot(xs, hsCompositeOuter, 'linewidth', 2, 'color', 'black');
%     plot(xsJet, hsJet, 'linewidth', 2);
%     hold off;
    plot(xs, hsCompositeOuter);
    hold on;
%     plot(xsJet, hsInnerLower, '-o');
%     plot(xsJet, hsJet, '-o');
%     plot(xsJet, hsOverlapJet, '-o');
    plot(xsJet, hsCompositeJet);
    hold off;


    xlim([epsilon * d - 0.005, epsilon * d + 0.015]);
    ylim([0, 0.02]);
    pbaspect([1 1 1]);
    
    drawnow;
    pause(0.1);
   
end
