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
set(0, 'defaultFigureRenderer', 'painters')
set(groot, 'DefaultLegendInterpreter', 'latex');



%% Load in color map
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;
blueCol = cmap(1, :);
redCol = cmap(end, :);

%% Load parameters
tmax = 1;
ts = linspace(0, tmax, 1e3)';

% Substrate parameters
[epsilon, L, q, omega] = plateparameters(); % Substrate parameters

% Stationary substrate coefficients
zeroTerm = zeros(size(ts));
StationarySubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm);

% Flat substrate coefficients
[ws, w_ts, w_tts] = flatsubstrate(ts, q, omega);
SubstrateCoefficients = substratecoefficients(ws, w_ts, w_tts);

%% Determine time dependents (apart from force)
% Stationary substrate time dependents
StationaryTimeDependents = timedependents(ts, StationarySubstrateCoefficients);

% Flat substrate time dependents
TimeDependents = timedependents(ts, SubstrateCoefficients);

%% Determine forces
% Stationary substrate forces
[Stationary_Fs_composite, Stationary_Fs_outer, Stationary_Fs_inner] ...
    = substrateforce(ts, StationaryTimeDependents, epsilon);

% Flat substrate time dependents
[Fs_composite, Fs_outer, Fs_inner] ...
    = substrateforce(ts, TimeDependents, epsilon);


%% Determine energies
% Stationary substrate energies
[Stationary_Es_outer, Stationary_Es_jets] ...
    = dropletenergy(ts, StationaryTimeDependents, StationarySubstrateCoefficients, epsilon);

% Stationary substrate energies
[Es_outer, Es_splash_sheet] ...
    = dropletenergy(ts, TimeDependents, SubstrateCoefficients, epsilon);

%% Substrate motion plot
% tiledlayout(1, 3);
% 
% % Substrate position
% nexttile;
% plot(ts, ws, 'color', 'black', 'linewidth', 2);
% xlabel("$t$");
% ylabel("$w(t$)");
% grid on;
% 
% % Substrate velocity
% nexttile;
% plot(ts, w_ts / epsilon, 'color', 'black', 'linewidth', 2);
% xlabel("$t$");
% ylabel("$w'(t) / \epsilon$");
% grid on;
% 
% % Substrate acceleration
% nexttile;
% plot(ts, w_tts / epsilon^2, 'color', 'black', 'linewidth', 2);
% xlabel("$t$");
% ylabel("$w''(t) / \epsilon^2$");
% grid on;
% 
% % Figure settings
% set(gcf,'position', [100, 100, 960, 300]);
% pause(0.1);
% 
% % Export figure
% exportgraphics(gcf,'png/AnalyticalImposed.png', 'Resolution', 300);
% exportgraphics(gcf,'eps/AnalyticalImposed.eps', 'Resolution', 300);

%% Turnover point plot
figNo = 2;
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryTimeDependents.ds, 'color', 'black', 'linewidth', 2);
plot(ts, TimeDependents.ds, 'color', redCol, 'linewidth', 2);

legend(["Stationary substrate", "Flat substrate"], ...
    'Location', 'southeast');
xlabel("$t$");
ylabel("$d_0(t)$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
exportgraphics(gca,'png/TurnoverPointsAxi.png', 'Resolution', 300);
exportgraphics(gca,'eps/TurnoverPointsAxi.eps', 'Resolution', 300);

%% Jet thickness
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryTimeDependents.Js, 'color', 'black', 'linewidth', 2);
plot(ts, TimeDependents.Js, 'color', redCol, 'linewidth', 2);

legend(["Stationary substrate", "Flat substrate"], ...
    'Location', 'northwest');
xlabel("$t$");
ylabel("$J(t)$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
exportgraphics(gca,'png/JetThicknessAxi.png', 'Resolution', 300);
exportgraphics(gca,'eps/JetThicknessAxi.eps', 'Resolution', 300);

%% Force plot

figure(figNo);
figNo = figNo + 1;
hold on;

% Plot composite forces
plot(ts, Stationary_Fs_composite, 'color', 'black', 'linewidth', 2);
plot(ts, Fs_composite, 'color', redCol, 'linewidth', 2);

% Plot outer forces
plot(ts, Stationary_Fs_outer, 'color', 'black', 'linewidth', 2, 'linestyle', '--');
plot(ts, Fs_outer, 'color', redCol, 'linewidth', 2, 'linestyle', '--');

lh = legend(["Stationary substrate", "Flat substrate"], ...
    'Location', 'northwest');

xlabel("$t$");
ylabel("$F(t)$");

ylim([0, 1.5]);
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
exportgraphics(gca,'png/ForceAxi.png', 'Resolution', 300);
exportgraphics(gca,'eps/ForceAxi.eps', 'Resolution', 300);

%% Energy plot

figure(figNo);
figNo = figNo + 1;
hold on;

% Plot outer energies
plot(ts, Stationary_Es_outer, 'color', 'black', 'linewidth', 1.5);
plot(ts, Es_outer, 'color', redCol, 'linewidth', 1.5);

% Plot jet energies
plot(ts, Stationary_Es_jets, 'color', 'black', 'linewidth', 3, 'linestyle', '--');
plot(ts, Es_splash_sheet, 'color', redCol, 'linewidth', 3, 'linestyle', '--');

legend(["Stationary substrate", "Flat substrate"], ...
    'Location', 'northwest');
xlabel("$t$");
ylabel("$E_{K}(t)$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
exportgraphics(gca,'png/EnergyAxi.png', 'Resolution', 300);
exportgraphics(gca,'eps/EnergyAxi.eps', 'Resolution', 300);


