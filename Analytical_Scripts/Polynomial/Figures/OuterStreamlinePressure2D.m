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

addpath("../");

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Load in color map
% mapObj = load("red_blue_cmap.mat");
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;

%% Parameters
epsilon = 0.1;
t = 1;

%% Substrate definition
% Oscillatory time dependence
L = 1;
q = 0.1;
omega = 2.5;

a = @(t) q * L^2 * (t^2 + 1 - cos(omega * t));
a_t = @(t) q * L^2 * (2 * t + omega * sin(omega * t));
a_tt = @(t) q * L^2 * (2 + omega^2 * cos(omega * t));

b = @(t) - a(t) / L^2;
b_t = @(t) - a_t(t) / L^2;
b_tt = @(t) - a_tt(t) / L^2;

%% Save structure arrays
% Stationary substrate case
StationarySubstrateCoefficients = substratecoefficients(0, 0, 0, 0, 0, 0, epsilon);
StationaryTimeDependents = timedependents(t, StationarySubstrateCoefficients);

% Quadratic substrate case
SubstrateCoefficients = substratecoefficients(a(t), b(t), a_t(t), ...
    b_t(t), a_tt(t), b_tt(t), epsilon);
TimeDependents = timedependents(t, SubstrateCoefficients);

%% Plotting parameters
d_stat = StationaryTimeDependents.ds;

% Spatial quantities
noPoints = 1e3; % Number of grid points per dimension
xMax = 2 * d_stat; % Maximum value of x
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
% grid on;

%% Creates mesh grid of values
[X, Z] = meshgrid(xs, zs);
Zeta = X + 1i * Z; % Complex variable on the mesh grid 

%% Pressure surf plot
% Determines right hand side coordinates
xsPos = linspace(0, xMax, noPoints / 2);
[XPos, ZPos] = meshgrid(xsPos, zs);
ZetaPos = XPos + 1i * ZPos;

% Determines stationary substrate solution in the right-hand-side
d_stat = StationaryTimeDependents.ds;
A_stat = StationaryTimeDependents.As;
[~, ps_stat] = outersolution(ZetaPos, d_stat, A_stat, StationarySubstrateCoefficients);

% Find min and max of pressure
pMin = min(min(ps_stat));
pMax = max(max(ps_stat));

% Find logarithmic scaling for levels
levels = exp(linspace(-6, 1.0, noFillConts));

% % Plot the scaled pressure on the right hand side
contourf(XPos, ZPos, ps_stat, levels, 'Edgecolor', 'None');

% Plot stationary pressure on left-hand-side
xsNeg = linspace(-xMax, 0, noPoints / 2);
[XNeg, ZNeg] = meshgrid(xsNeg, zs);
ZetaNeg = XNeg + 1i * ZNeg;
PNeg = flip(ps_stat, 2);
contourf(XNeg, ZNeg, PNeg, levels, 'Edgecolor', 'None');

% Determines moving substrate solution on right-hand-side
d = TimeDependents.ds;
A = TimeDependents.As;
[Ws, ps] = outersolution(ZetaPos, d, A, SubstrateCoefficients);

% Plot the moving substrate pressure on the right hand side
contourf(XPos, ZPos, ps, levels, 'Edgecolor', 'None')

%% Plot streamline contours
% Stationary substrate
[W_stat, ~] = outersolution(ZetaNeg, d_stat, A_stat, StationarySubstrateCoefficients);

% Calculate streamfunction
Psi_stat = imag(W_stat);

% Plot contours of streamfunction (i.e. the streamlines)
contour(XNeg, ZNeg, Psi_stat, 15, 'Color', 0.25 * [1 1 1], 'Linewidth', 2)

% Moving substrate
[W, ~] = outersolution(ZetaPos, d, A, SubstrateCoefficients);

% Calculate streamfunction
Psi = imag(W);

% Plot contours of streamfunction (i.e. the streamlines)
contour(XPos, ZPos, Psi, 15, 'Color', 0.25 * [1 1 1], 'Linewidth', 2)


%% Scatter for turnover points
scatter(-d_stat, 0, [], 'black', 'filled');
scatter(d, 0, [], 'black', 'filled');

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
% set(gca, 'xtick',[-1, 0, 1]);
% xNames = {'$-d_0(t)$'; '0'; '$d_0(t)$'};
% set(gca, 'XTickLabel', xNames);

% y-axis settings
set(gca, 'ytick',[0, d_stat, 2 * d_stat]);
% yNames = {'0'; '$d_0(t)$'; '$2 d_0(t)$'};
% set(gca, 'YTickLabel', yNames);

box on;




% Axes line width
% axes1.LineWidth = 1.5;

% Set figure size
set(gcf,'position', [0, 0, 1200, 600]);

% Set rendered to Painters (incredibly slow but makes the figures better)
% set(gcf, 'Renderer', 'Painters');

%% Create figures
% Export png
% exportgraphics(gca,'png/OuterStreamlinePressure2D.png', 'Resolution', 300);

% Export eps (doesn't look great so left commented out)
% exportgraphics(gca,'eps/OuterStreamlinePressure2D.eps', 'Resolution', 300);