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
[epsilon, L, q, omega] = quadraticparameters(); 

tmax = 1;
ts = tmax;
xs = linspace(-L, L, 1e3);

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

%% Free-surface plot

% Array of SubstrateCoefficient structs
SubstrateCoefficientsArray = [StationarySubstrateCoefficients, ...
       FlatSubstrateCoefficients, QuadSubstrateCoefficients];

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
    SubstrateCoefficients = SubstrateCoefficientsArray(typeIdx);
    
    % Determine time dependents
    TimeDependents = timedependents(ts, SubstrateCoefficients);
    
    % Load turnover point
    d = TimeDependents.ds;
    
    % Create xs
    xs_Free_Surface = linspace(epsilon * d, L, 1e3);
    
    % Determine free-surface
    hs = outerfreesurface(xs_Free_Surface, ts, d, SubstrateCoefficients, epsilon);
    
    % Plot free-surface
    h(typeIdx) = plot(xs_Free_Surface, hs, 'linewidth', 2, 'color', lineColor, ...
        'Displayname', type);
    plot(-xs_Free_Surface, hs, 'linewidth', 2, 'color', lineColor, ...
        'Displayname', type);
    
    %% Plot the substrate
    a = SubstrateCoefficients.aHats;
    b = SubstrateCoefficients.bHats / epsilon^2;
    plot(xs, -epsilon^2 * (a * ones(size(xs)) + b * xs.^2), ...
        'Color', lineColor, 'Linestyle', ':', 'Linewidth', 2);
end


%% Figure settings
xlim([-L, L]);
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
