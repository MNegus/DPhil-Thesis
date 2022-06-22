%%CompositePressure2D.m
%   Script to make a figure comparing the composite pressure to the outer
%   and inner pressures in 2D.

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
[epsilon, k, q, omega] = substrateparameters(); 

xs = linspace(0, 0.4, 5e3);

t = 2;

%% Load in substrate function
SubstrateFunctions = substratefunctions("stationary");

%% Find pressures
[ps_composite, ps_outer, ps_inner] ...
    = substratepressure(xs, t, SubstrateFunctions, epsilon);

% Restrict outer solution where zero
ps_outer(ps_outer == 0) = nan;

%% Set line colour
lineColor = 'Black';

%% Plot pressures
figure(1);
hold on;

% Composite pressure
h(3) = plot(xs, ps_composite, 'Linewidth', 2, 'Color', lineColor, ...
        'Displayname', 'Composite solution');

% Outer pressure
h(1) = plot(xs, ps_outer, 'Linewidth', 3, 'Color', redCol, ...
    'Linestyle', '--', 'Displayname', 'Outer solution');

% Inner pressure
h(2) = plot(xs, ps_inner, 'Linewidth', 3, 'Color', blueCol, ...
    'Linestyle', ':', 'Displayname', 'Inner solution');


%% Plot turnover line
xline(epsilon * SubstrateFunctions.d(t), 'Linewidth', 0.5, 'Color', 0.5 * [1 1 1], ...
    'Linestyle', '--');

%% Figure options
grid on;
box on;

% Axes labels
xlabel("$x$");
ylabel("$p(x, 0, t)$");

% Set axes limits
xlim([0, 0.35]);
ylim([0, 40]);
% xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25]);

% Legend
legend(h(1:3), 'Location', 'Northwest');

% Set figure position
set(gcf,'position', [100, 100, 600, 400]);

set(gcf, 'Renderer', 'Painters');
pause(0.1);

filename = "CompositePressure2D";
savefig(gcf, sprintf("fig/%s.fig", filename));
exportgraphics(gca, sprintf("png/%s.png", filename), 'Resolution', 300);
exportgraphics(gca,sprintf("eps/%s.eps", filename), 'Resolution', 300);

