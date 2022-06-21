%% OuterPressureAxi
% Script to produce a figure of the streamlines and pressure in the outer
% region for the axisymmetric impact case. 

clear;
close all;

addpath("../");

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Load in color map
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;
colormap(cmap);

%% Parameter definitions
noFillConts = 50;

% Substrate parameters
[epsilon, L, q, omega] = plateparameters(); 

tmax = 0.25;
ts = tmax;
centre_rs = linspace(0, L, 1e4);

% Stationary substrate coefficients
zeroTerm = zeros(size(ts));
StationarySubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm);
StationaryTimeDependents = timedependents(ts, StationarySubstrateCoefficients);

% Flat substrate coefficients
[ws, w_ts, w_tts] = flatsubstrate(ts, q, omega);
SubstrateCoefficients = substratecoefficients(ws, w_ts, w_tts);
TimeDependents = timedependents(ts, SubstrateCoefficients);

%% Save structure arrays
A_stationary = StationaryTimeDependents.As;

A_moving = TimeDependents.As;

%% Plotting parameters
minP = 0;
maxP = 5000;
levels = linspace(minP, maxP, 100);

%% Surf plot for pressure
numRs = 1e2;
numThetas = 1e3;
centre_rs = sqrt(sqrt(linspace(0, 1, numRs)));
thetas = linspace(0, 2 * pi, numThetas);

figure(1);
hold on;

% Stationary substrate pressure (left)
thetaLeft = linspace(pi, 2 * pi, numThetas / 2);
[centreR, Theta] = meshgrid(centre_rs, thetaLeft);
R = centreR .* sin(Theta);
Z = 1 - centreR .* cos(Theta);
Pleft = Pfun(A_stationary, R, Z);
Pleft(Pleft < 0) = 0;
contourf(R, Z, Pleft, levels, 'Edgecolor', 'None');

% Moving substrate pressure (right)
thetaRight = linspace(0, pi, numThetas / 2);
[centreR, Theta] = meshgrid(centre_rs, thetaRight);
R = centreR .* sin(Theta);
Z = 1 - centreR .* cos(Theta);
Pright = Pfun(A_moving, R, Z);
Pright(Pright < 0) = 0;
contourf(R, Z, Pright, levels, 'Edgecolor', 'None');

plot(sin(thetas), 1 + cos(thetas), 'Color', 'Black');

xlim([-1, 1]);
ylim([0, 2]);
pbaspect([1 1 1]);

% Colourbar 
cb = colorbar('Location', 'Northoutside');
cb.Label.String = '$P_0(r, z, t)$';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
cb.XTick = [0, 2500, 5000];
caxis([minP maxP]);

%% Whole figure settings
grid on;
xlabel("$r$");
ylabel("$z$");
set(gcf,'position', [0, 0, 800, 600]);

% Set position of colour bar to be above the axes at the centre (a bit
% hacky but works for now)
axesPos = get(gca,'position');
figPos = get(gcf, 'Position');
cbPos = get(cb,'Position');

cbPos(3) = 0.75 * cbPos(3); % Half the width of the colour bar
cbPos(1) = (1.3 * axesPos(3) / 8) + 0.5 * cbPos(3); % Set left of the colour bar to be 1/8 along the axis
cbPos(2) = 0.83;

set(cb,'Position',cbPos)
set(gca,'Position',axesPos);

% Axes limits
pertb = 0.1;
extent = 1 + pertb;

% Add white dividing line
xline(0, 'color', 'white');

% Axes limits
xlim([-extent, extent]);
ylim([0, 2 * extent]);

% Set position of colour bar to be above the axes at the centre (a bit
% hacky but works for now)
axesPos = get(gca,'position');
figPos = get(gcf, 'Position');
cbPos = get(cb,'Position');

cbPos(3) = 0.75 * cbPos(3); % Half the width of the colour bar
cbPos(1) = (2.25 * axesPos(3) / 8) + 0.5 * cbPos(3); % Set left of the colour bar to be 1/8 along the axis
cbPos(2) = 0.83;

set(cb,'Position',cbPos)
set(gca,'Position',axesPos);

% Export figure
% set(gcf, 'Renderer', 'Painters');
% exportgraphics(gca,'png/OuterOuterWhole2D.png', 'Resolution', 300);



%% Wedge plot
plotLim = 0.1;
rMin = 1 - 2 * plotLim;
ang = asin(plotLim / rMin);
numRs = 5e2;
numThetas = 5e2;

centre_rs = sqrt(linspace((1 - 2 * plotLim).^2, 1, numRs));
thetas = linspace(-ang, ang, numThetas);

levels = linspace(minP, maxP, 30);

figure(2);
hold on;
colormap(cmap);

% Stationary substrate pressure (left)
thetaLeft = linspace(-ang, 0, numThetas / 2);
[centreR, Theta] = meshgrid(centre_rs, thetaLeft);
R = centreR .* sin(Theta);
Z = 1 - centreR .* cos(Theta);
Pleft = Pfun(A_stationary, R, Z);
Pleft(Pleft < 0) = 0;
contourf(R, Z, Pleft, levels, 'Edgecolor', 'None');

% Moving substrate pressure (right)
thetaRight = linspace(0, ang, numThetas / 2);
[centreR, Theta] = meshgrid(centre_rs, thetaRight);
R = centreR .* sin(Theta);
Z = 1 - centreR .* cos(Theta);
Pright = Pfun(A_moving, R, Z);
Pright(Pright < 0) = 0;
contourf(R, Z, Pright, levels, 'Edgecolor', 'None');

% Add black line for interface
plot(sin(thetas), 1 - cos(thetas), 'Color', 'Black');

% Add white dividing line
xline(0, 'color', 'white');

xlim([-plotLim, plotLim]);
ylim([0, 2 * plotLim]);
pbaspect([1 1 1]);
grid on;
xlabel("$r$");
ylabel("$z$");

set(gcf,'position', [0, 0, 300, 300]);
% set(gcf, 'Renderer', 'Painters');
% exportgraphics(gca,'png/OuterOuterZoomed2D.png', 'Resolution', 300);


function Ps = Pfun(A, R, Z)
    Ps = A * (1 - R.^2 - (Z - 1).^2) ./ (2 * (R.^2 + Z.^2).^(3/2));
end
