%% OuterPressureAxi
% Script to produce a figure of the streamlines and pressure in the outer
% region for the 2D impact case. Appears in Section 3.3.6, under the Wagner
% theory chapter. 

clear;
close all;

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
noFillConts = 100;

% Time dependent quantities
d = 1; % Let the turnover point be at 1
t = d^2 / 4; % Solves for t as a function of d
d_t = 1 / sqrt(t); % Time derivative of d
G = 3 * sqrt(3) * t; % Factor in front of pressure

%% Surf plot for pressure
numRs = 5e2;
numThetas = 1e3;
radials = sqrt(sqrt(linspace(0, 1, numRs)));
thetas = linspace(0 , 2 * pi, 1e3);

[radials, Theta] = meshgrid(radials, thetas);
R = radials .* sin(Theta);
Z = 1 - radials .* cos(Theta);

P = PFun(R, Z, G);
P(P < 0) = 0;


minP = 0;
maxP = 1000;
max(max(P))
% levels = exp(linspace(-10, log(maxP), noFillConts));
levels = linspace(minP, maxP, 100);
% levels = 30;


figure(1);
hold on;
contourf(R, Z, P, levels, 'Edgecolor', 'none');
% surf(X, Z, P, 'Edgecolor', 'None');
% surf(X, Z, P)
% % contourf(X, Z, log(P), levels);
plot(sin(thetas), 1 + cos(thetas), 'Color', 'Black');
xlim([-1, 1]);
ylim([0, 2]);
pbaspect([1 1 1]);

% Colourbar 
cb = colorbar('Location', 'Northoutside');
cb.Label.String = '$P_0(x, z, t)$';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
cb.XTick = [0, 50, 100, 150, 200, 250];
caxis([minP maxP]);

%% Whole figure settings
grid on;
xlabel("$x$");
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
xlim([-extent, extent]);
ylim([0, 2 * extent]);

% Export figure
% set(gcf, 'Renderer', 'Painters');
% exportgraphics(gca,'png/OuterOuterWhole2D.png', 'Resolution', 300);



%% Wedge plot
plotLim = 0.1;
rMin = 1 - 2 * plotLim;
ang = asin(plotLim / rMin);
% ang = pi / 10;
thetas = linspace(-ang, ang, 5e2);
radials = sqrt(linspace((1 - 2 * plotLim).^2, 1, 5e2));

[radials, Theta] = meshgrid(radials, thetas);
R = radials .* sin(Theta);
Z = 1 - radials .* cos(Theta);

P = PFun(R, Z, G);
P(P < 0) = 0;

% levels = exp(linspace(-10, log(maxP), 5 * noFillConts));
levels = linspace(minP, maxP, 30);

close(figure(2));
figure(2);
colormap(cmap);
hold on;
contourf(R, Z, P, levels, 'Edgecolor', 'black');
plot(sin(thetas), 1 - cos(thetas), 'Color', 'Black');

xlim([-plotLim, plotLim]);
ylim([0, 2 * plotLim]);
pbaspect([1 1 1]);
grid on;
xlabel("$x$");
ylabel("$z$");

set(gcf,'position', [0, 0, 300, 300]);
% set(gcf, 'Renderer', 'Painters');
% exportgraphics(gca,'png/OuterOuterZoomed2D.png', 'Resolution', 300);

function Ps = PFun(R, Z, G)
    Ps = G * (1 - R.^2  - (Z - 1).^2) ./ (2 * (R.^2 + Z.^2).^(3/2));
end