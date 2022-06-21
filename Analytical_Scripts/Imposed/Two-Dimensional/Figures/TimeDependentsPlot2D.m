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
[epsilon, L, q, omega] = quadraticparameters(); % Substrate parameters

% Stationary substrate coefficients
zeroTerm = zeros(size(ts));
StationarySubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm, zeroTerm, zeroTerm, zeroTerm, epsilon);

% Flat substrate coefficients
[aFlats, a_tFlats, a_ttFlats, bFlats, b_tFlats, b_ttFlats] ...
    = flatsubstrate(ts, q, omega); 
FlatSubstrateCoefficients ...
    = substratecoefficients(aFlats, bFlats, a_tFlats, b_tFlats, a_ttFlats, b_ttFlats, epsilon);


% Quadratic substrate coefficients
[aQuads, a_tQuads, a_ttQuads, bQuads, b_tQuads, b_ttQuads] ...
    = quadraticsubstrate(ts, L, q, omega);
QuadSubstrateCoefficients ...
    = substratecoefficients(aQuads, bQuads, a_tQuads, b_tQuads, a_ttQuads, b_ttQuads, epsilon);

%% Determine time dependents (apart from force)
% Stationary substrate time dependents
StationaryTimeDependents = timedependents(ts, StationarySubstrateCoefficients);

% Flat substrate time dependents
FlatTimeDependents = timedependents(ts, FlatSubstrateCoefficients);

% Quadratic substrate time dependents
QuadTimeDependents = timedependents(ts, QuadSubstrateCoefficients);

%% Determine forces
% Stationary substrate forces
[Stationary_Fs_composite, Stationary_Fs_outer, Stationary_Fs_inner] ...
    = substrateforce(ts, StationaryTimeDependents, epsilon);

% Flat substrate time dependents
[Flat_Fs_composite, Flat_Fs_outer, Flat_Fs_inner] ...
    = substrateforce(ts, FlatTimeDependents, epsilon);

% Quadratic substrate time dependents
[Quad_Fs_composite, Quad_Fs_outer, Quad_Fs_inner] ...
    = substrateforce(ts, QuadTimeDependents, epsilon);

%% Determine energies
% Stationary substrate energies
[Stationary_Es_outer, Stationary_Es_jets] ...
    = dropletenergy(ts, StationaryTimeDependents, StationarySubstrateCoefficients, epsilon);

% Stationary substrate energies
[Flat_Es_outer, Flat_Es_jets] ...
    = dropletenergy(ts, FlatTimeDependents, FlatSubstrateCoefficients, epsilon);

% Stationary substrate energies
[Quad_Es_outer, Quad_Es_jets] ...
    = dropletenergy(ts, QuadTimeDependents, QuadSubstrateCoefficients, epsilon);

%% Substrate motion plot
tiledlayout(1, 3);

% Substrate position
nexttile;
ws = aFlats;
plot(ts, ws, 'color', 'black', 'linewidth', 2);
xlabel("$t$");
ylabel("$w(t$)");
grid on;

% Substrate velocity
nexttile;
w_ts = a_tFlats;
plot(ts, w_ts / epsilon, 'color', 'black', 'linewidth', 2);
xlabel("$t$");
ylabel("$w'(t) / \epsilon$");
grid on;

% Substrate acceleration
nexttile;
w_tts = a_ttFlats;
plot(ts, w_tts / epsilon^2, 'color', 'black', 'linewidth', 2);
xlabel("$t$");
ylabel("$w''(t) / \epsilon^2$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 960, 300]);
pause(0.1);

% Export figure
exportgraphics(gcf,'png/AnalyticalImposed.png', 'Resolution', 300);
exportgraphics(gcf,'eps/AnalyticalImposed.eps', 'Resolution', 300);

%% Turnover point plot
figNo = 2;
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryTimeDependents.ds, 'color', 'black', 'linewidth', 2);
plot(ts, FlatTimeDependents.ds, 'color', redCol, 'linewidth', 2);
plot(ts, QuadTimeDependents.ds, 'color', blueCol, 'linewidth', 2);

legend(["Stationary substrate", "Flat substrate", "Curved substrate"], ...
    'Location', 'southeast');
xlabel("$t$");
ylabel("$d_0(t)$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
exportgraphics(gca,'png/TurnoverPoints2D.png', 'Resolution', 300);
exportgraphics(gca,'eps/TurnoverPoints2D.eps', 'Resolution', 300);

%% Jet thickness
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryTimeDependents.Js, 'color', 'black', 'linewidth', 2);
plot(ts, FlatTimeDependents.Js, 'color', redCol, 'linewidth', 2);
plot(ts, QuadTimeDependents.Js, 'color', blueCol, 'linewidth', 2);

legend(["Stationary substrate", "Flat substrate", "Curved substrate"], ...
    'Location', 'northwest');
xlabel("$t$");
ylabel("$J(t)$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
exportgraphics(gca,'png/JetThickness2D.png', 'Resolution', 300);
exportgraphics(gca,'eps/JetThickness2D.eps', 'Resolution', 300);

%% Force plot

figure(figNo);
figNo = figNo + 1;
hold on;

% Plot composite forces
plot(ts, Stationary_Fs_composite, 'color', 'black', 'linewidth', 2);
plot(ts, Flat_Fs_composite, 'color', redCol, 'linewidth', 2);
plot(ts, Quad_Fs_composite, 'color', blueCol, 'linewidth', 2);

% Plot outer forces
plot(ts, Stationary_Fs_outer, 'color', 'black', 'linewidth', 2, 'linestyle', '--');
plot(ts, Flat_Fs_outer, 'color', redCol, 'linewidth', 2, 'linestyle', '--');
plot(ts, Quad_Fs_outer, 'color', blueCol, 'linewidth', 2, 'linestyle', '--');

lh = legend(["Stationary substrate", "Flat substrate", "Curved substrate"], ...
    'Location', 'best');
lh.Position(1) = 0.5 - lh.Position(3)/2; 
lh.Position(2) = 0.575 - lh.Position(4)/2;

xlabel("$t$");
ylabel("$F(t)$");

ylim([0, 8]);
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
exportgraphics(gca,'png/Force2D.png', 'Resolution', 300);
exportgraphics(gca,'eps/Force2D.eps', 'Resolution', 300);

%% Energy plot

figure(figNo);
figNo = figNo + 1;
hold on;

% Plot outer energies
plot(ts, Stationary_Es_outer, 'color', 'black', 'linewidth', 1.5);
plot(ts, Flat_Es_outer, 'color', redCol, 'linewidth', 1.5);
plot(ts, Quad_Es_outer, 'color', blueCol, 'linewidth', 1.5);

% Plot jet energies
plot(ts, Stationary_Es_jets, 'color', 'black', 'linewidth', 3, 'linestyle', '--');
plot(ts, Flat_Es_jets, 'color', redCol, 'linewidth', 3, 'linestyle', '--');
plot(ts, Quad_Es_jets, 'color', blueCol, 'linewidth', 3, 'linestyle', '--');

legend(["Stationary substrate", "Flat substrate", "Curved substrate"], ...
    'Location', 'northwest');
xlabel("$t$");
ylabel("$E_{K}(t)$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 600, 300]);
pause(0.1);

% Export figure
exportgraphics(gca,'png/Energy2D.png', 'Resolution', 300);
exportgraphics(gca,'eps/Energy2D.eps', 'Resolution', 300);


