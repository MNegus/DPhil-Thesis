%%PressureEvolution2D.m
%   Script to make figures comparing pressure (composite) along the
%   substrate in time for the stationary, plate and quadratic substrate
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
ts = linspace(0.05, tmax, 1000)';
xMax = 2 * epsilon * sqrt(tmax);
xs = linspace(0, xMax, 1e4);

%% Load in substrate functions
StationaryFunctions = substratefunctions("stationary");
FlatFunctions = substratefunctions("flat");
CurvedFunctions = substratefunctions("curved");


%% Pressure in time plot

% Array of TimeDependents structs
FunctionsArray = [StationaryFunctions, FlatFunctions, CurvedFunctions];

% Array of types
typeArr = ["stationary", "flat", "curved"];

freq = 100; % How frequent the lines are

tiledlayout(1, 3);

tileNo = 1;
for typeIdx = 1 : length(typeArr)
    type = typeArr(typeIdx);
    
    SubstrateFunctions = FunctionsArray(typeIdx);
    
    % Move to next tile
    nexttile(tileNo);
    hold on;
    tileNo = tileNo + 1;
    
    for k = 1 : freq : length(ts)
        t = ts(k);
        if type == "stationary"
            lineColor = 'black';
            
        elseif type == "flat"
            lineColor = redCol;
            
        else
            lineColor = blueCol;
        end

        % Load in composite pressure
        [ps, ~, ~] = substratepressure(xs, t, SubstrateFunctions, epsilon);

        % Plot pressure
        plot(xs, ps, 'Linewidth', 1, 'color', lineColor);
        
    end
    
    %% Plot maximum pressure
    [pMaxs, xMaxs] = pressuremax(ts, SubstrateFunctions, epsilon);
    plot(xMaxs, pMaxs, 'color', lineColor, 'linestyle', '--');
    
    %% Tile figure settings
    xlim([0, xMax]);
    ylim([-10, 300]);
    xlabel("$x$");
    ylabel("$p(x, t)$");
    grid on;
    
end

%% Figure settings
set(gcf,'position', [100, 100, 960, 300]);
pause(0.1);

% set(gcf, 'Renderer', 'Painters');

% Export figure
% savefig(gcf, 'fig/PressureEvolution2D.fig');
% exportgraphics(gcf,'png/PressureEvolution2D.png', 'Resolution', 300);
% exportgraphics(gcf,'eps/PressureEvolution2D.eps', 'Resolution', 300);