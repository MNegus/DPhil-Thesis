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

%% Load parameters
% Substrate parameters
[epsilon, k, q, omega] = substrateparameters(); 

tmax = 1;
t = tmax;

L = epsilon * 2 * sqrt(tmax);
xs = linspace(-1.1 * L, 1.1 * L, 1e3);

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
typeArr = ["Stationary substrate solution", ...
    "Flat substrate solution", "Quadratic substrate solution"];

figure(1);
hold on;
for typeIdx = 1 : length(typeArr)
    type = typeArr(typeIdx);
    
    if type == "Stationary substrate solution"
       lineColor = 'black';
    elseif type == "Flat substrate solution"
        lineColor = redCol;
    else
        lineColor = blueCol;
    end
    
    % Save SubstrateCoefficients
    SubstrateFunctions = SubstrateFunctionsArray(typeIdx);
    
    % Load turnover point
    d = SubstrateFunctions.d(t);
    
    % Create xs
    xs_Free_Surface = linspace(epsilon * d, 2 * L, 1e3);
    
    % Determine free-surface
    hHats = outerfreesurface(xs_Free_Surface, t, SubstrateFunctions);
    hs = epsilon^2 * hHats;
    
    % Plot free-surface
    h(typeIdx) = plot(xs_Free_Surface, hs, 'linewidth', 2, 'color', lineColor, ...
        'Displayname', type);
    plot(-xs_Free_Surface, hs, 'linewidth', 2, 'color', lineColor, ...
        'Displayname', type);
    
    %% Plot the substrate
    ws = SubstrateFunctions.w(xs, t);
    plot(xs, -epsilon^2 * ws, ...
        'Color', lineColor, 'Linestyle', ':', 'Linewidth', 2);
end


%% Figure settings
% xlim([-L, L]);
legend(h(1:3), 'Location', 'North');

grid on;

xlabel("$x$");
ylabel("$z$");

set(gcf,'position', [100, 100, 700, 350]);
set(gcf, 'Renderer', 'painters');

% Export figure
savefig(gcf, 'fig/FreeSurface2D.fig');
exportgraphics(gca,'png/FreeSurface2D.png', 'Resolution', 300);
exportgraphics(gca,'eps/FreeSurface2D.eps', 'Resolution', 300);
