close all;
clear;
addpath("../FreeSurface");

SubstrateFunctions = imposedsubstratefunctions("stationary", "axi");

%% Parameters
ts = linspace(0.001, 10, 1e2);

tiledlayout(1, 2);
for t = ts
% for t = 0.5
    
    d = SubstrateFunctions.d(t);
    epsilon = SubstrateFunctions.epsilon;
    xMaxUpper = 1;
    xMaxLower = 1.5 * epsilon * d;

    %% Outer outer solution
    [xsLower, hsLower, xsUpper, hsUpper] = outerouterfreesurface(t, SubstrateFunctions);
    plot(xsLower, hsLower);
    hold on;
    plot(xsUpper, hsUpper);
    hold off;
    
    xlim([0, 2]);
    ylim([0, 2]);
    pbaspect([1 1 1]);
    %% Load composite solution
%     [xs, hsTurnover, hsFull] = outer_inner_jet_freesurface_composite(xMaxUpper, xMaxLower, t, SubstrateFunctions);
    
    %% Plot solutions
%     nexttile(1);
%     plot(xs, hsFull);
%     xlim([0, 1]);
%     ylim([0, 1]);  
%     pbaspect([1 1 1]);
%     
%     nexttile(2);
%     plot(xs, hsTurnover);
%     xlim([epsilon * d - 0.005, epsilon * d + 0.015]);
%     ylim([0, 0.02]);
%     pbaspect([1 1 1]);
%     
    drawnow;
    pause(0.1);
   
end
