%%FreeSurface2D.m
%   Plots the free surface for the stationary, flat and quadratic substrate
%   cases.

clear;
close all;

addpath("../");
addpath("../Pressures");

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Load in color map
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;
blueCol = cmap(1, :);
redCol = cmap(end, :);
blackCol = [0, 0, 0];

%% Load parameters
% Substrate parameters
[epsilon, k, q, omega] = substrateparameters(); 

tmax = 1;
% ts = linspace(1e-10, tmax, 10);
ts = [0.5, tmax];


L = epsilon * 2 * sqrt(tmax);
% xs = linspace(-1.1 * L, 1.1 * L, 1e3);
xs = linspace(0, 2 * L, 1e3);

%% Load in substrate functions
dimension = "2D";
StationaryFunctions = substratefunctions("stationary", dimension);
FlatFunctions = substratefunctions("flat", dimension);
CurvedFunctions = substratefunctions("curved", dimension);

%% Free-surface plot

% Array of SubstrateCoefficient structs
SubstrateFunctionsArray = [StationaryFunctions, ...
       FlatFunctions, CurvedFunctions];

% Array of types
types = ["stationary", "flat", "curved"];
colors = [blackCol; redCol; blueCol];


figNo = 1;
for t = ts
    figure(figNo);
    figNo = figNo + 1;
    for typeIdx = 1 : length(types)
        type = types(typeIdx);

        % Save SubstrateCoefficients
        SubstrateFunctions = substratefunctions(types(typeIdx), dimension);
        lineColor = colors(typeIdx, :);

        % Load turnover point
        d = SubstrateFunctions.d(t);

        % Create xs
        xHats_Free_Surface = linspace(d, 2 * L / epsilon, 1e3);

        % Determine free-surface
        hHats = outerfreesurface(xHats_Free_Surface, t, SubstrateFunctions);
        hs = hHats;

        % Plot free-surface
        h(typeIdx) = plot(xHats_Free_Surface, hs, 'linewidth', 2, 'color', lineColor, ...
            'Displayname', type);
        hold on;

        %% Plot the substrate
        ws = SubstrateFunctions.w(xs, t);
        plot(xs / epsilon, -ws, ...
            'Color', lineColor, 'Linestyle', ':', 'Linewidth', 2);
    end
    hold off;
    
    ylim([-2, 6]);
    
    legend(h(1:length(types)), 'Location', 'North');

    grid on;

    xlabel("$x$");
    ylabel("$z$");

    set(gcf,'position', [100, 100, 500, 400]);
end



%% Figure settings
% xlim([-L, L]);

set(gcf, 'Renderer', 'painters');

% Export figure
% savefig(gcf, 'fig/FreeSurface2D.fig');
% exportgraphics(gca,'png/FreeSurface2D.png', 'Resolution', 300);
% exportgraphics(gca,'eps/FreeSurface2D.eps', 'Resolution', 300);
