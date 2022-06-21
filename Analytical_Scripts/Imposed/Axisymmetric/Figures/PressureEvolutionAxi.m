%%PressureEvolutionAxi.m
%   Script to make figures comparing pressure (composite) along the
%   substrate in time for the stationary and plate cases.
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
[epsilon, L, q, omega] = plateparameters(); % Substrate parameters

tmax = 1;
ts = linspace(1e-6, tmax, 1e3)';
rs = linspace(0, 0.2, 1e4);

% Stationary substrate coefficients
zeroTerm = zeros(size(ts));
StationarySubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm);
StationaryTimeDependents = timedependents(ts, StationarySubstrateCoefficients);

% Flat substrate coefficients
[ws, w_ts, w_tts] = flatsubstrate(ts, q, omega);
SubstrateCoefficients = substratecoefficients(ws, w_ts, w_tts);
TimeDependents = timedependents(ts, SubstrateCoefficients);

%% Pressure in time plot

% Array of TimeDependents structs
TimeDependentsArray = [StationaryTimeDependents, TimeDependents];

% Array of types
typeArr = ["stationary", "flat"];

freq = 100; % How frequent the lines are

tiledlayout(1, 2);

tileNo = 1;
for typeIdx = 1 : length(typeArr)
    type = typeArr(typeIdx)
    
    % Move to next tile
    nexttile(tileNo);
    hold on;
    tileNo = tileNo + 1;
    
    for k = freq : freq : length(ts)
        t = ts(k);
        if type == "stationary"
            SubstrateCoefficients ...
                = substratecoefficients(0, 0, 0);
            lineColor = 'black';
            
        elseif type == "flat"
            SubstrateCoefficients = substratecoefficients(ws(k), w_ts(k), w_tts(k));
            lineColor = redCol;
        end
        
        % Load in time dependents
        TimeDependents = timedependents(t, SubstrateCoefficients);

        % Load in composite pressure
        [ps, ~, ~] ...
            = substratepressure(rs, TimeDependents, epsilon);

        % Plot pressure
        plot(rs, ps, 'Linewidth', 1, 'color', lineColor);
        
    end
    
    %% Plot maximum pressure
    [pMaxs, xMaxs] = pressuremax(TimeDependentsArray(typeIdx), epsilon);
    plot(xMaxs, pMaxs, 'color', lineColor, 'linestyle', '--');
    
    %% Tile figure settings
    xlim([0, 0.2]);
    ylim([-10, 150]);
    xlabel("$r$");
    ylabel("$p(r, t)$");
    grid on;
    
end

%% Figure settings
set(gcf,'position', [100, 100, 640, 300]);
pause(0.1);

% set(gcf, 'Renderer', 'Painters');
% 
% % Export figure
% exportgraphics(gcf,'png/PressureEvolution2D.png', 'Resolution', 300);
% exportgraphics(gcf,'eps/PressureEvolution2D.eps', 'Resolution', 300);