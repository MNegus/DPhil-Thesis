%% OuterStreamlinePressure
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
mapObj = load("red_blue_cmap.mat");
cmap = mapObj.cmap;
colormap(cmap);

%% Parameter definitions
noFillConts = 100;

% Time dependent quantities
d = 1; % Let the turnover point be at 1
t = d^2 / 4; % Solves for t as a function of d
d_t = 1 / sqrt(t); % Time derivative of d
A = d * d_t; % Factor in front of pressure
D = d^2 / 2;

%% Surf plot for pressure
rs = linspace(0, 1, 1e3);
thetas = linspace(0, 2 * pi, 1e3);
[R, Theta] = meshgrid(rs, thetas);
X = R .* cos(Theta);
Z = 1 + R .* sin(Theta);

P = - A * (X.^2 + (Z - 1).^2 - 1) ./ (2 * (X.^2 + Z.^2));

% log(min(min(P)))
log(max(max(P)))
% Find logarithmic scaling for levels
levels = exp(linspace(-6.5, 4.5, noFillConts));

figure(1);
hold on;
% contourf(X, Z, P, levels, 'Edgecolor', 'None');
plot(cos(thetas), 1 + sin(thetas), 'Color', 'Black');
xlim([-1, 1]);
ylim([0, 2]);
pbaspect([1 1 1]);

%% Streamlines
Psi = X ./ (X.^2 + Z.^2);
contour(X, Z, Psi, 15, 'Color', 0.25 * [1 1 1], 'Linewidth', 2)

%% Colourbar 
% cb = colorbar('Location', 'Northoutside');