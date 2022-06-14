%% JetStreamlinePressure2D.m
% Script to produce a figure of the streamlines and pressure in the jet
% region for the 2D impact case. Appears in Section 3.3.5, under the Wagner
% theory chapter. 
% 
% For the figure, we neglect any motion of the membrane, taking w == 0. 
%

clear;
close all;

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

mapObj = load("red_blue_cmap.mat");
cmap = mapObj.cmap;

figure(1);
colormap(cmap);
hold on;

% Axes labels
xlabel('$\bar{x}$');
ylabel('$\bar{z}$');

box on;


%% Parameters
epsilon = 1;

% Time
t = 1 / 4;

% Time dependent functions
d = @(t) 2 * sqrt(t);
d_t = @(t) 1 ./ sqrt(t);
d_tt = @(t) -1 ./ (2 * t.^(3/2));
d_ttt = @(t) 3 ./ (4 * t.^(5/2));
J = @(t) pi * d(t) ./ (8 * d_t(t).^2);

% Solution for x
x = @(tau) 2 * d_t(tau) .* (t - tau) + d(tau);

% Solution for u and derivatives
u = @(tau) 2 * d_t(tau);
dudx = @(tau) 2 * d_tt(tau) ./ (2 * d_tt(tau) .* (t - tau) - d_t(tau));
detJ = @(tau) d_t(tau) - 2 * d_tt(tau) .* (t - tau);
d2udx2 = @(tau) 2 * (d_t(tau) .* d_tt(tau) - 3 * d_tt(tau)) ./ detJ(tau).^2;
d2udxdt = @(tau) 2 * (d_tt(tau).^2 - d_t(tau) .* d_ttt(tau)) ./ detJ(tau).^2;

% Solution for h
h = @(tau) (d_t(tau) .* J(tau)) ./ (d_t(tau) - 2 * d_tt(tau) .* (t - tau));

% Solution for p
p = @(tau, z) epsilon^4 * (u(tau) .* d2udx2(tau) .* (h(tau) - z) ...
    + 0.5 * (d2udxdt(tau) - dudx(tau).^2) .* (h(tau).^2 - z.^2));


% Range for x
xMin = d(t);
xMax = 3 * d(t);
zMax = 1.5 * h(t);

%% Finds range for tau
% Minimum tau is where x = xMax
tauMin = fsolve(@(tau) xMax - x(tau), 1e-6);
taus = linspace(tauMin, t, 1e3);
xs = x(taus);

%% Plots the pressure
zScales = linspace(0, 1, 1e3);

% Create a scaled meshgrid
[Taus, ZScales] = meshgrid(taus, zScales);

Hs = h(Taus);
Xs = x(Taus);
Zs = Hs .* ZScales;
Ps = p(Taus, Zs);
contourf(Xs, Zs, Ps, 10, 'Edgecolor', 'None');

%% Plot the free surface
hs = h(taus);

figure(1);
plot(xs, hs, 'color', 'black', 'linewidth', 2);

%% Create quiver plot of velocity
us = u(taus);
% xQuiv = xMin : 0.1 : xMax;
tausQuiv = linspace(tauMin, t, 5);
zQuiv = 0 : 0.01 : h(t);
[TausQuiv, ZQuiv] = meshgrid(tausQuiv, zQuiv);
XQuiv = x(TausQuiv);
UQuiv = u(TausQuiv);
VQuiv = 0 * u(TausQuiv);

% quiver(XQuiv, ZQuiv, UQuiv, VQuiv);

%% Plot colour bar for the pressure plot
cb = colorbar('Location', 'Northoutside');
cb.Label.String = '$\bar{p}_4(\bar{x}, \bar{z}, t)$';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';

% Set position of colour bar to be above the axes at the centre (a bit
% hacky but works for now)
x1 = get(gca,'position');
x = get(cb,'Position');
x(3) = 0.5 * x(3);
x(2) = 0.99 * x(2);
x(1) = 2.25 * x(1);
set(cb,'Position',x)
set(gca,'position',x1)



%% Figure properties
xlim([xMin, xMax]);
% ylim([0, zMax]);
% pbaspect([1 zMax / (xMax- xMin) 1]);

% x-axis settings
set(gca, 'xtick',[1, 2, 3]);
xNames = {'$d_0(t)$'; '$2 d_0(t)$'; '$3 d_0(t)$'};
set(gca, 'XTickLabel', xNames);

% y-axis settings
set(gca, 'ytick',[0, 0.5 * J(t), J(t)]);
yNames = {'$0$'; '$0.5 \, J(t)$'; '$J(t)$'};
set(gca, 'YTickLabel', yNames);

% Set figure size
set(gcf,'position', [0, 0, 800, 400]);


%% Create figures
% Export png
exportgraphics(gca,'png/JetStreamlinePressure2D.png', 'Resolution', 300);