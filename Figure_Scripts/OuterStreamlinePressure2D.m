%% OuterStreamlinePressure
% Script to produce a figure of the streamlines and pressure in the outer
% region for the 2D impact case. Appears in Section 3.3.3, under the Wagner
% theory chapter. 
% 
% For the figure, we neglect any motion of the membrane, taking w == 0. In
% this case, the complex potential is given by
% W(\zeta) =  \phi + i \psi = i(\zeta^2 - d(t)^2)^{1/2}, 
% and the pressure is
% p = \Re(i A(t) / (\zeta^2 - d(t)^2)^{1/2}), A(t) = d(t) d'(t).
%
% We plot the streamlines as a contour plot for constant \psi, and then a
% colour plot for the pressure p. We pick a generic value of d(t) for
% plotting clarity. 

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

%% Parameter definitions
% Time dependent quantities
d = 1; % Let the turnover point be at 1
t = d^2 / 4; % Solves for t as a function of d
d_t = 1 / sqrt(t); % Time derivative of d
A = d * d_t; % Factor in front of pressure

% Spatial quantities
noPoints = 1e3; % Number of grid points per dimension
xMax = 2 * d; % Maximum value of x
zMax = xMax; % Maximum value of z

xs = linspace(-xMax, xMax, noPoints); % Array of points in x
zs = linspace(0, zMax, noPoints); % Array of points in z

% Filled contour options
noFillConts = 30; % Number of filled contours
scaleFun = @(P) log(P); % Scale function for the filled contours
inverseScaleFun = @(val) exp(val); % Inverse of scale function

%% Figure properties
figure1 = figure();
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Set colour map to cool-warm
colormap(cmap);

% Axes labels
xlabel('$\hat{x}$');
ylabel('$\hat{z}$');

% Show axes grid
grid on;

%% Creates mesh grid of values
[X, Z] = meshgrid(xs, zs);
Zeta = X + 1i * Z; % Complex variable on the mesh grid 

%% Pressure surf plot
% Determine the pressure in the right-hand side
xsPos = linspace(0, xMax, noPoints / 2);
[XPos, ZPos] = meshgrid(xsPos, zs);
ZetaPos = XPos + 1i * ZPos;
PPos = real(1i * A ./ (ZetaPos.^2 - d^2).^0.5);

% Find min and max of pressure
pMin = min(min(PPos));
pMax = max(max(PPos));

% Find logarithmic scaling for levels
levels = exp(linspace(-6, 1.8, noFillConts));

% Plot the scaled pressure on the right hand side
contourf(XPos, ZPos, PPos, levels, 'Edgecolor', 'None');

% Left-hand figure, made by reflecting the right-hand figure about the line
% x = 0
xsNeg = linspace(-xMax, 0, noPoints / 2);
[XNeg, ZNeg] = meshgrid(xsNeg, zs);
ZetaNeg = XNeg + 1i * ZNeg;
PNeg = flip(PPos, 2);
contourf(XNeg, ZNeg, PNeg, levels, 'Edgecolor', 'None');

%% Plot streamline contours
% Calculate streamfunction
Psi = imag(1i * (Zeta.^2 - d^2).^0.5);

% Plot contours of streamfunction (i.e. the streamlines)
contour(X, Z, Psi, 15, 'Color', 'black', 'Linewidth', 2)

%% Plot colour bar for the pressure plot
cb = colorbar('Location', 'Northoutside');
cb.Label.String = '$\hat{p}_0(\hat{x}, \hat{z}, t)$';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';

% Set position of colour bar to be above the axes at the centre (a bit
% hacky but works for now)
x1 = get(gca,'position');
x = get(cb,'Position');
x(3) = 0.5 * x(3);
x(2) = 0.95 * x(2);
x(1) = 2.475 * x(1);
set(cb,'Position',x)
set(gca,'position',x1)


%% Figure settings
pbaspect([1 zMax / (2 * xMax) 1]); % Aspect ratio of plot
set(gca, 'Layer', 'Top'); % Set axes to be on top layer

% x-axis settings
set(gca, 'xtick',[-1, 0, 1]);
xNames = {'$-d_0(t)$'; '0'; '$d_0(t)$'};
set(gca, 'XTickLabel', xNames);

% y-axis settings
set(gca, 'ytick',[0, 1, 2]);
yNames = {'0'; '$d_0(t)$'; '$2 d_0(t)$'};
set(gca, 'YTickLabel', yNames);

% Axes line width
axes1.LineWidth = 1.5;

% Set figure size
set(gcf,'position', [0, 0, 1200, 600]);

% Set rendered to Painters (incredibly slow but makes the figures better)
set(gcf, 'Renderer', 'Painters');

%% Create figures
% Export png
exportgraphics(gca,'png/OuterStreamlinePressure2D.png', 'Resolution', 300);

% Export eps (doesn't look great so left commented out)
% exportgraphics(gca,'eps/OuterStreamlinePressure2D.eps', 'Resolution', 300);