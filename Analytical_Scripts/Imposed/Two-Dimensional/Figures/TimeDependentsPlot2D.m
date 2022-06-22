%%TimeDependents2D.m
%   Script to make figures comparing the time-dependent quantities (i.e.
%   substrate position, turnover point, jet thickness etc. for a stationary
%   substrate, a flat moving and a quadratic moving substrate.

clear;
close all;

addpath("../");
addpath("../Forces");
addpath("../Energies");

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Load in color map
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;
blueCol = cmap(1, :);
redCol = cmap(end, :);

%% Load parameters

% Substrate parameters
[epsilon, L, q, omega] = substrateparameters(); % Substrate parameters

tmax = 1;
ts = linspace(0, tmax, 1e3)';

%% Load in substrate functions
StationaryFunctions = substratefunctions("stationary");
FlatFunctions = substratefunctions("flat");
CurvedFunctions = substratefunctions("curved");

%% Determine forces
% Stationary substrate forces
[Stationary_Fs_composite, Stationary_Fs_outer, Stationary_Fs_inner] ...
    = substrateforce(ts, StationaryFunctions, epsilon);

% Flat substrate time dependents
[Flat_Fs_composite, Flat_Fs_outer, Flat_Fs_inner] ...
    = substrateforce(ts, FlatFunctions, epsilon);

% Quadratic substrate time dependents
[Curved_Fs_composite, Curved_Fs_outer, Curved_Fs_inner] ...
    = substrateforce(ts, CurvedFunctions, epsilon);

%% Determine energies
% Stationary substrate energies
[Stationary_Es_outer, Stationary_Es_jets] ...
    = dropletenergy(ts, StationaryFunctions, epsilon);

% Stationary substrate energies
[Flat_Es_outer, Flat_Es_jets] ...
    = dropletenergy(ts, FlatFunctions, epsilon);

% Stationary substrate energies
[Curved_Es_outer, Curved_Es_jets] ...
    = dropletenergy(ts, CurvedFunctions, epsilon);

%% Substrate motion plot
fig = tiledlayout(1, 3);

% Substrate position
nexttile;
ws = FlatFunctions.a(ts);
plot(ts, ws, 'color', 'black', 'linewidth', 2);
xlabel("$t$");
ylabel("$w(t$)");
grid on;

% Substrate velocity
nexttile;
w_ts = FlatFunctions.a_t(ts);
plot(ts, w_ts / epsilon, 'color', 'black', 'linewidth', 2);
xlabel("$t$");
ylabel("$w'(t) / \epsilon$");
grid on;

% Substrate acceleration
nexttile;
w_tts = FlatFunctions.a_tt(ts);
plot(ts, w_tts / epsilon^2, 'color', 'black', 'linewidth', 2);
xlabel("$t$");
ylabel("$w''(t) / \epsilon^2$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 960, 300]);
pause(0.1);

% Export figure
filename = "AnalyticalImposed";
savefig(gcf, sprintf("fig/%s.fig", filename));
exportgraphics(gcf, sprintf("png/%s.png", filename), 'Resolution', 300);
exportgraphics(gcf,sprintf("eps/%s.eps", filename), 'Resolution', 300);

%% Turnover point plot
figNo = 2;
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryFunctions.d(ts), 'color', 'black', 'linewidth', 2);
plot(ts, FlatFunctions.d(ts), 'color', redCol, 'linewidth', 2);
plot(ts, CurvedFunctions.d(ts), 'color', blueCol, 'linewidth', 2);

legend(["Stationary substrate", "Flat substrate", "Curved substrate"], ...
    'Location', 'southeast');
xlabel("$t$");
ylabel("$d_0(t)$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
filename = "TurnoverPoints2D";
savefig(gcf, sprintf("fig/%s.fig", filename));
exportgraphics(gca, sprintf("png/%s.png", filename), 'Resolution', 300);
exportgraphics(gca,sprintf("eps/%s.eps", filename), 'Resolution', 300);

%% Jet thickness
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryFunctions.J(ts), 'color', 'black', 'linewidth', 2);
plot(ts, FlatFunctions.J(ts), 'color', redCol, 'linewidth', 2);
plot(ts, CurvedFunctions.J(ts), 'color', blueCol, 'linewidth', 2);

legend(["Stationary substrate", "Flat substrate", "Curved substrate"], ...
    'Location', 'northwest');
xlabel("$t$");
ylabel("$J(t)$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
filename = "JetThickness2D";
savefig(gcf, sprintf("fig/%s.fig", filename));
exportgraphics(gca, sprintf("png/%s.png", filename), 'Resolution', 300);
exportgraphics(gca,sprintf("eps/%s.eps", filename), 'Resolution', 300);

%% Force plot

figure(figNo);
figNo = figNo + 1;
hold on;

% Plot composite forces
plot(ts, Stationary_Fs_composite, 'color', 'black', 'linewidth', 2);
plot(ts, Flat_Fs_composite, 'color', redCol, 'linewidth', 2);
plot(ts, Curved_Fs_composite, 'color', blueCol, 'linewidth', 2);

% Plot outer forces
plot(ts, Stationary_Fs_outer, 'color', 'black', 'linewidth', 2, 'linestyle', '--');
plot(ts, Flat_Fs_outer, 'color', redCol, 'linewidth', 2, 'linestyle', '--');
plot(ts, Curved_Fs_outer, 'color', blueCol, 'linewidth', 2, 'linestyle', '--');

lh = legend(["Stationary substrate", "Flat substrate", "Curved substrate"], ...
    'Location', 'northwest');
% lh.Position(1) = 0.52 - lh.Position(3)/2; 
% lh.Position(2) = 0.575 - lh.Position(4)/2;

xlabel("$t$");
ylabel("$F(t)$");

ylim([0, 15]);
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
filename = "Force2D";
savefig(gcf, sprintf("fig/%s.fig", filename));
exportgraphics(gca, sprintf("png/%s.png", filename), 'Resolution', 300);
exportgraphics(gca,sprintf("eps/%s.eps", filename), 'Resolution', 300);

%% Energy plot

figure(figNo);
figNo = figNo + 1;
hold on;

% Plot outer energies
plot(ts, Stationary_Es_outer, 'color', 'black', 'linewidth', 1.5);
plot(ts, Flat_Es_outer, 'color', redCol, 'linewidth', 1.5);
plot(ts, Curved_Es_outer, 'color', blueCol, 'linewidth', 1.5);

% Plot jet energies
plot(ts, Stationary_Es_jets, 'color', 'black', 'linewidth', 3, 'linestyle', '--');
plot(ts, Flat_Es_jets, 'color', redCol, 'linewidth', 3, 'linestyle', '--');
plot(ts, Curved_Es_jets, 'color', blueCol, 'linewidth', 3, 'linestyle', '--');

legend(["Stationary substrate", "Flat substrate", "Curved substrate"], ...
    'Location', 'northwest');
xlabel("$t$");
ylabel("$E_{K}(t)$");
grid on;
ylim([0, 0.035]);

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
filename = "Energy2D";
savefig(gcf, sprintf("fig/%s.fig", filename));
exportgraphics(gca, sprintf("png/%s.png", filename), 'Resolution', 300);
exportgraphics(gca,sprintf("eps/%s.eps", filename), 'Resolution', 300);

