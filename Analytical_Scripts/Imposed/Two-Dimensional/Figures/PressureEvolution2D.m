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
[epsilon, L, q, omega] = quadraticparameters(); 

tmax = 1;
ts = linspace(0.05, tmax, 1000)';
xs = linspace(0, L, 1e4);

% Stationary substrate coefficients
zeroTerm = zeros(size(ts));
StationarySubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm, zeroTerm, zeroTerm, zeroTerm, epsilon);

% Flat substrate coefficients
[aFlats, a_tFlats, a_ttFlats, bFlats, b_tFlats, b_ttFlats] ...
    = flatsubstrate(ts, q, omega); 
FlatSubstrateCoefficients ...
    = substratecoefficients(aFlats, bFlats, a_tFlats, b_tFlats, a_ttFlats, b_ttFlats, epsilon);


% Quadratic substrate coefficients
[aQuads, a_tQuads, a_ttQuads, bQuads, b_tQuads, b_ttQuads] ...
    = quadraticsubstrate(ts, L, q, omega);
QuadSubstrateCoefficients ...
    = substratecoefficients(aQuads, bQuads, a_tQuads, b_tQuads, a_ttQuads, b_ttQuads, epsilon);

%% Determine time dependents
% Stationary substrate time dependents
StationaryTimeDependents = timedependents(ts, StationarySubstrateCoefficients);

% Flat substrate time dependents
FlatTimeDependents = timedependents(ts, FlatSubstrateCoefficients);

% Quadratic substrate time dependents
QuadTimeDependents = timedependents(ts, QuadSubstrateCoefficients);


%% Pressure in time plot

% Array of TimeDependents structs
TimeDependentsArray = [StationaryTimeDependents, FlatTimeDependents, QuadTimeDependents];

% Array of types
typeArr = ["stationary", "flat", "quad"];

freq = 100; % How frequent the lines are

tiledlayout(1, 3);

tileNo = 1;
for typeIdx = 1 : length(typeArr)
    type = typeArr(typeIdx)
    
    % Move to next tile
    nexttile(tileNo);
    hold on;
    tileNo = tileNo + 1;
    
    for k = 1 : freq : length(ts)
        t = ts(k);
        if type == "stationary"
            SubstrateCoefficients ...
                = substratecoefficients(0, 0, 0, 0, 0, 0, epsilon);
            lineColor = 'black';
            
        elseif type == "flat"
            SubstrateCoefficients ...
                = substratecoefficients(aFlats(k), bFlats(k), ...
                a_tFlats(k), b_tFlats(k), a_ttFlats(k), b_ttFlats(k), ...
                epsilon);
            lineColor = redCol;
            
        else
            SubstrateCoefficients ...
                = substratecoefficients(aQuads(k), bQuads(k), ...
                a_tQuads(k), b_tQuads(k), a_ttQuads(k), b_ttQuads(k), ...
                epsilon);
            lineColor = blueCol;
        end
        
        
        % Load in time dependents
        TimeDependents = timedependents(t, SubstrateCoefficients);

        % Load in composite pressure
        [ps, ~, ~] ...
            = substratepressure(xs, SubstrateCoefficients, TimeDependents, epsilon);

        % Plot pressure
        plot(xs, ps, 'Linewidth', 1, 'color', lineColor);
        
    end
    
    %% Plot maximum pressure
    [pMaxs, xMaxs] = pressuremax(TimeDependentsArray(typeIdx), epsilon);
    plot(xMaxs, pMaxs, 'color', lineColor, 'linestyle', '--');
    
    %% Tile figure settings
    xlim([0, 0.2]);
    ylim([-10, 150]);
    xlabel("$x$");
    ylabel("$p(x, t)$");
    grid on;
    
end

%% Figure settings
set(gcf,'position', [100, 100, 960, 300]);
pause(0.1);

set(gcf, 'Renderer', 'Painters');

% Export figure
savefig(gcf, 'fig/PressureEvolution2D.fig');
exportgraphics(gcf,'png/PressureEvolution2D.png', 'Resolution', 300);
exportgraphics(gcf,'eps/PressureEvolution2D.eps', 'Resolution', 300);