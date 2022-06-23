%%PressureEvolution2D.m
%   Script to make figures comparing pressure (composite) along the
%   substrate in time for the stationary, plate and quadratic substrate
%   cases.

clear;
close all;

addpath("../");
addpath("../Pressures");

painters = true;

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

if painters
    set(0, 'defaultFigureRenderer', 'painters');
else
    set(0, 'defaultFigureRenderer', 'opengl');
end

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
noPoints = 1e3;

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

tiledlayout(3, 1);

tileNo = 1;
for typeIdx = 1 : length(typeArr)
    type = typeArr(typeIdx);
    
    % Set line colors
    if type == "stationary"
            lineColor = 'black';
            
    elseif type == "flat"
        lineColor = redCol;

    else
        lineColor = blueCol;
    end
    
    % Load substrate functions
    SubstrateFunctions = FunctionsArray(typeIdx);
    
    % Move to next tile
    nexttile(tileNo);
    hold on;
    tileNo = tileNo + 1;
    
    for k = 1 : freq : length(ts)
        t = ts(k);

        %% Determine x values, clustered around the maximum pressure point
        % Find maximum pressure point
        [~, x_pMax] = pressuremax(t, SubstrateFunctions, epsilon);
        
        fineWidth = 1e-3;
        leftWidth = x_pMax/ xMax;
        rightWidth = (xMax - x_pMax) / xMax;
        xsLower = linspace(0, x_pMax - fineWidth, floor(leftWidth * noPoints));
        xsFine = linspace(x_pMax - fineWidth, x_pMax + fineWidth, 100);
            
        xsUpper = linspace(x_pMax + fineWidth, xMax, floor(rightWidth * noPoints));
        
        xs = [xsLower, xsFine, xsUpper];
        
        % Load in composite pressure
        [ps, ~, ~] = substratepressure(xs, t, SubstrateFunctions, epsilon);

        % Plot pressure
        plot(xs, ps, 'Linewidth', 1, 'color', lineColor);
        
    end
    
    %% Plot maximum pressure
    tLongs = linspace(0, 2 * tmax);
    [pMaxs, xMaxs] = pressuremax(tLongs, SubstrateFunctions, epsilon);
    plot(xMaxs, pMaxs, 'color', lineColor, 'linestyle', '--');
    
    %% Tile figure settings
    xlim([0, xMax]);
    ylim([-10, 300]);
    xlabel("$x$");
%     ylabel("$p_{\textrm{comp}}(x, t)$");
    ylabel("$p_{{comp}}(x, t)$");
    grid on;
    
end

%% Figure settings
set(gcf,'position', [100, 100, 800, 600]);
pause(0.1);

% Export figure
filename = "PressureEvolution2D";
savefig(gcf, sprintf("fig/%s.fig", filename));
exportgraphics(gcf, sprintf("png/%s.png", filename), 'Resolution', 300);
exportgraphics(gcf,sprintf("eps/%s.eps", filename), 'Resolution', 300);
