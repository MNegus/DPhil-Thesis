%%FreeSurfaceAxi.m
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
[epsilon, L, q, omega] = quadraticparameters(); 

tmax = 1;
ts = tmax;
rs = linspace(0, L, 1e3);

% Stationary substrate coefficients
zeroTerm = zeros(size(ts));
StationarySubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm);

% Stationary substrate coefficients
zeroTerm = zeros(size(ts));
StationarySubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm);

% Flat substrate coefficients
[ws, w_ts, w_tts] = flatsubstrate(ts, q, omega);
SubstrateCoefficients = substratecoefficients(ws, w_ts, w_tts);


%% Free-surface plot

% Array of SubstrateCoefficient structs
SubstrateCoefficientsArray = [StationarySubstrateCoefficients, ...
       SubstrateCoefficients];

% Array of types
typeArr = ["Stationary substrate solution", ...
    "Flat substrate solution"];

figure(1);
hold on;
for typeIdx = 1 : length(typeArr)
    type = typeArr(typeIdx);
    
    if type == "Stationary substrate solution"
       lineColor = 'black';
    elseif type == "Flat substrate solution"
        lineColor = redCol;
    end
    
    % Save SubstrateCoefficients
    SubstrateCoefficients = SubstrateCoefficientsArray(typeIdx);
    
    % Determine time dependents
    TimeDependents = timedependents(ts, SubstrateCoefficients);
    
    % Load turnover point
    d = TimeDependents.ds;
    
    % Create xs
    rs_Free_Surface = linspace(epsilon * d, L, 1e3);
    
    % Determine free-surface
    hs = outerfreesurface(rs_Free_Surface, ts, d, SubstrateCoefficients, epsilon);
    
    % Plot free-surface
    h(typeIdx) = plot(rs_Free_Surface, hs, 'linewidth', 2, 'color', lineColor, ...
        'Displayname', type);
    
    %% Plot the substrate
    w = SubstrateCoefficients.ws;
    plot(rs, -epsilon^2 * w * ones(size(rs)), ...
        'Color', lineColor, 'Linestyle', ':', 'Linewidth', 2);
end


%% Figure settings
xlim([0, L]);
legend(h(1:2), 'Location', 'Northwest');

grid on;

xlabel("$r$");
ylabel("$z$");

set(gcf,'position', [100, 100, 600, 350]);
set(gcf, 'Renderer', 'painters');

% Export figure
savefig(gcf, 'fig/FreeSurfaceAxi.fig');
exportgraphics(gcf,'png/FreeSurfaceAxi.png', 'Resolution', 300);
exportgraphics(gcf,'eps/FreeSurfaceAxi.eps', 'Resolution', 300);
