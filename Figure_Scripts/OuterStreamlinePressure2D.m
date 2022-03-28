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
xMax = 1.5 * d; % Maximum value of x
zMax = xMax; % Maximum value of z

xs = linspace(-xMax, xMax, noPoints); % Array of points in x
zs = linspace(0, zMax, noPoints); % Array of points in z

%% Figure properties
figure(1);
hold on;

%% Creates mesh grid of values
[X, Z] = meshgrid(xs, zs);
Zeta = X + 1i * Z;

%% Pressure surf plot
noConts = 20;
scaleFun = @(P) 2 * log(P);

xsPos = linspace(0, xMax, noPoints / 2);
[XPos, ZPos] = meshgrid(xsPos, zs);
ZetaPos = XPos + 1i * ZPos;
PPos = real(1i * A ./ (ZetaPos.^2 - d^2).^0.5);

contourf(XPos, ZPos, scaleFun(PPos), noConts, 'Edgecolor', 'None')

xsNeg = linspace(-xMax, 0, noPoints / 2);
[XNeg, ZNeg] = meshgrid(xsNeg, zs);
ZetaNeg = XNeg + 1i * ZNeg;
PNeg = flip(PPos, 2);

contourf(XNeg, ZNeg, scaleFun(PNeg), noConts, 'Edgecolor', 'None')

colormap(cmap);

%% Plot streamline contours
Psi = imag(1i * (Zeta.^2 - d^2).^0.5);
contour(X, Z, Psi, 15, 'Color', 'black')

%% Figure settings
pbaspect([1 zMax / (2 * xMax) 1]);
%%
exportgraphics(gca,'OuterStreamlinePressure2D.png', 'Resolution', 300)