%%CompositeForce2D.m
%   Script to make a figure comparing the composite force to the
%   leading-order force in 2D.

clear;
close all;

addpath("../Forces");

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Parameters
epsilon = 0.1;
dmax = 0.2; % Maximum turnover point value
% tmax = (dmax / (2 * epsilon))^2
tmax = 1
ts = linspace(0, tmax, 1e3);

%% Load in plate positions and time-dependents
ws = zeros(size(ts));
w_ts = zeros(size(ts));
w_tts = zeros(size(ts));

PlatePositions = platepositions(ws, w_ts, w_tts);
TimeDependents = timedependents(ts, PlatePositions);

%% Determine forces
[Fs_composite, Fs_outer, Fs_inner] ...
    = substrateforce(ts, PlatePositions, TimeDependents, epsilon);
Fs_composite(1)
Fs_outer(1)

%% Plot forces
figure(1);
hold on;
plot(ts, Fs_outer, 'Linewidth', 3, 'Color', 'black', ...
        'Linestyle', '--', 'Displayname', 'Leading-order');
plot(ts, Fs_composite, 'Linewidth', 1.85, 'Color', 'Black', ...
        'Displayname', 'Composite');

%% Figure options
grid on;
box on;

% Axes labels
xlabel("$t$");
ylabel("$F(t)$");

% Set axes limits
% xlim([0, 0.26]);
ylim([5.6, 6.6]);

% Legend
legend('Location', 'north');

% Set figure position
set(gcf,'position', [100, 100, 960, 400]);

set(gcf, 'Renderer', 'Painters');

% Export
exportgraphics(gca,'png/CompositeForce2D.png', 'Resolution', 300);
exportgraphics(gca,'eps/CompositeForce2D.eps', 'Resolution', 300);


