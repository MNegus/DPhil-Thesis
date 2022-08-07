%% PlateDisplacementComparison.m
%

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Forces/");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;

redCol = cmap(end, :);
blueCol = cmap(1, :);

%% Figure options
fontsize = 11;
lineWidth = 1.25;
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize', fontsize);
set(0, 'DefaultTextFontSize', fontsize);
set(0,'defaultLegendFontSize', fontsize, 'DefaultLegendFontSizeMode','manual');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(0, 'DefaultLegendInterpreter', 'latex');

%% Parameters
% Computational parameters
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
MAX_TIMESTEP = T_MAX / DELTA_T;

% Master DNS dir
dns_master_dir = "/home/michael/scratch/DPhil_DNS_Data/";

% Analytical parameters
epsilon = 1;
tMax = 1.5 * (1 / 3); % A bit after the stationary max value

% Plate parameters
ALPHA = 2;
BETA = 0;
GAMMA = 20;

%% Load stationary plate DNS solutions
stat_dir = sprintf("%s/Stationary_Plate/axi", dns_master_dir);
stat_output_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", stat_dir));
ts = stat_output_mat(:, 1) - IMPACT_TIME;
FsStat = stat_output_mat(:, 3);

%% Load moving plate DNS solutions
param_dir = sprintf("%s/Moving_Plate/ALPHA-%g_BETA-%g_GAMMA-%g", dns_master_dir, ...
    ALPHA, BETA, GAMMA);
moving_output_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", param_dir));
ts = moving_output_mat(:, 1) - IMPACT_TIME;
FsMoving = moving_output_mat(:, 2);
wsMoving = moving_output_mat(:, 6);
w_tsMoving = moving_output_mat(:, 7);
w_ttsMoving = moving_output_mat(:, 8);

%% Load outer solution
% Solve plate equation
[tsOuter, wsOuter, w_tsOuter, w_ttsOuter] ...
    = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "outer");

% Find substrate functions
OuterSubstrateFunctions = platesubstratefunctions(tsOuter, ...
    wsOuter, w_tsOuter, w_ttsOuter, epsilon);

% Find where turnover point reaches 1
dsOuter = OuterSubstrateFunctions.d(tsOuter);
tIdxMaxOuter = sum(dsOuter <= 1);

% Restrict solutions temporally
tsOuter = tsOuter(1 : tIdxMaxOuter);
wsOuter = wsOuter(1 : tIdxMaxOuter);
w_tsOuter = w_tsOuter(1 : tIdxMaxOuter);
w_ttsOuter = w_ttsOuter(1 : tIdxMaxOuter);
dsOuter = dsOuter(1 : tIdxMaxOuter);

% Force solution
FsOuter = outerforce(tsOuter, OuterSubstrateFunctions);

%% Load composite solution
% Solve plate equation
[tsComp, wsComp, w_tsComp, w_ttsComp] ...
    = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

% Find substrate functions
CompSubstrateFunctions = platesubstratefunctions(tsComp, ...
    wsComp, w_tsComp, w_ttsComp, epsilon);

% Find where turnover point reaches 1
dsComp = CompSubstrateFunctions.d(tsComp);
tIdxMaxComp = sum(dsComp <= 1);

% Restrict solutions temporally
tsComp = tsComp(1 : tIdxMaxComp);
wsComp = wsComp(1 : tIdxMaxComp);
w_tsComp = w_tsComp(1 : tIdxMaxComp);
w_ttsComp = w_ttsComp(1 : tIdxMaxComp);
dsComp = dsComp(1 : tIdxMaxComp);

% Force solution
[FsComp, ~, ~] ...
    = substrateforce(tsComp, CompSubstrateFunctions);

%% Load stationary solutions
tMax = 1 / 3; % Max such that 1 = d(t)
tsStat = linspace(0, tMax, 1e3);

% Find substrate functions
StatSubstrateFunctions = platesubstratefunctions(tsStat, ...
    zeros(size(tsStat)), zeros(size(tsStat)), zeros(size(tsStat)), epsilon);

% Find ds
dsStat = StatSubstrateFunctions.d(tsStat);

% Force solution
[FsCompStat, FsOuterStat, ~] ...
    = substrateforce(tsStat, StatSubstrateFunctions);

%% Plot substrate displacement
tiledlayout(3, 1);
width = 6;
height = 5;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

xTicks = [-0.125 : 0.125 : 0.75];

% Position
nexttile;
hold on;
xline(tsOuter(end), '--', 'color', 0.75 * [1 1 1]);
xline(tsComp(end), ':', 'color', 0.75 * [1 1 1]);
h(1) = plot(ts, -wsMoving, 'color', blueCol, 'linewidth', lineWidth);
h(2) = plot(tsOuter, -wsOuter, 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
h(3) = plot(tsComp, -wsComp, 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);

grid on;
box on;
xlabel("$t$");
xticks(xTicks);

ylabel("$-w(t)$");
ylim([-0.15, 0.02]);

% Velocity
nexttile;
hold on;
xline(tsOuter(end), '--', 'color', 0.75 * [1 1 1]);
xline(tsComp(end), ':', 'color', 0.75 * [1 1 1]);
plot(ts, -w_tsMoving, 'color', blueCol, 'linewidth', lineWidth);
plot(tsOuter, -w_tsOuter, 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
plot(tsComp, -w_tsComp, 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);

grid on;
box on;

xlabel("$t$");
xticks(xTicks);

ylabel("$-\dot{w}(t)$");
ylim([-0.4, 0.1]);

% Acceleration
nexttile;
hold on;
xline(tsOuter(end), '--', 'color', 0.75 * [1 1 1]);
xline(tsComp(end), ':', 'color', 0.75 * [1 1 1]);
plot(ts, -w_ttsMoving, 'color', blueCol, 'linewidth', lineWidth);
plot(tsOuter, -w_ttsOuter, 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
plot(tsComp, -w_ttsComp, 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);

grid on;
box on;

xlabel("$t$");
xticks(xTicks);

ylabel("$-\ddot{w}(t)$");
ylim([-1.5, 0.5]);

% Set figure options
lh = legend(h(1 : 3), ...
    ["DNS", "Analytical (leading-order)", "Analytical (composite)"], ...
    'NumColumns', 3);
lh.Layout.Tile = 'North'; 



set(gcf, 'Renderer', 'Painters');
pause(0.5);

figname = "PlateFigures/PlateComparison";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
saveas(gcf, sprintf("%s.png", figname));
% savefig(gcf, sprintf("%s.fig", figname));

%% Plot force
tiledlayout(2, 1);
width = 6;
height = 4;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Stationary solutions
nexttile;
hold on;
xline(tMax, '--', 'color', 0.75 * [1 1 1]);
xline(tMax, ':', 'color', 0.75 * [1 1 1]);
h(1) = plot(ts, FsStat, 'color', blueCol, 'linewidth', lineWidth);
h(2) = plot(tsStat, FsOuterStat, 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
h(3) = plot(tsStat, FsCompStat, 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);

grid on;
box on;

xlabel("$t$");
xticks(xTicks);

ylabel("$F(t)$");
ylim([0, 6]);

title("(a) Stationary substrate.", 'Fontsize', fontsize);
set(gca, 'TitleFontSizeMultiplier', 1);

% Moving solutions
nexttile;
hold on;
xline(tsOuter(end), '--', 'color', 0.75 * [1 1 1]);
xline(tsComp(end), ':', 'color', 0.75 * [1 1 1]);
plot(ts, FsMoving, 'color', blueCol, 'linewidth', lineWidth);
plot(tsOuter, FsOuter, 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
plot(tsComp, FsComp, 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);

grid on;
box on;

xlabel("$t$");
xticks(xTicks);

ylabel("$F(t)$");
ylim([0, 6]);

title("(b) Moving substrate.", 'Fontsize', fontsize);
set(gca, 'TitleFontSizeMultiplier', 1);


% Set figure options
lh = legend(h(1 : 3), ...
    ["DNS", "Analytical (leading-order)", "Analytical (composite)"], ...
    'NumColumns', 3);
lh.Layout.Tile = 'North'; 

set(gcf, 'Renderer', 'Painters');
pause(0.5);

figname = "PlateFigures/ForceComparison";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
saveas(gcf, sprintf("%s.png", figname));