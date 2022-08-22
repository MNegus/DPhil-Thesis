%% PressureStressComparison.m

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Pressures/");
addpath("../Analytical_Scripts/Forces/");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;

redCol = cmap(end, :);
blueCol = cmap(1, :);

%% Figure options
fontsize = 10;
lineWidth = 1.25;
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize', fontsize);
set(0, 'DefaultTextFontSize', fontsize);
set(0,'defaultLegendFontSize', fontsize, 'DefaultLegendFontSizeMode','manual');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(0, 'DefaultLegendInterpreter', 'latex');

%% Parameters
% Plate parameters
ALPHA = 2;
BETA = 0;
GAMMA = 20;
L = 2; 

% Computational parameters
DELTA_T = 1e-3;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
MAX_TIMESTEP = T_MAX / DELTA_T;

% DNS directories
dns_master_dir = "/home/michael/scratch/DPhil_DNS_Data/";
stat_dir = sprintf("%s/Stationary_Plate/axi", dns_master_dir);
moving_dir = sprintf("%s/Moving_Plate/ALPHA-%g_BETA-%g_GAMMA-%g", dns_master_dir, ...
    ALPHA, BETA, GAMMA);

% Analytical parameters
epsilon = 1;

% Analytical xs array
xsAnalytical = linspace(0, L, 1e3);

%% Load stationary analytical solutions
tMax = 1 / 3; % Max such that 1 = d(t)
tsStat = linspace(0, tMax, 1e3);

% Find substrate functions
StatSubstrateFunctions = platesubstratefunctions(tsStat, ...
    zeros(size(tsStat)), zeros(size(tsStat)), zeros(size(tsStat)), epsilon);

% Find ds
dsStat = StatSubstrateFunctions.d(tsStat);

%% Load moving analytical solutions
tMax = 1.5 * (1 / 3); % A bit after the stationary max value

% Solve plate equation
[tsComp, wsComp, w_tsComp, w_ttsComp] ...
    = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

% Find substrate functions
MovingSubstrateFunctions = platesubstratefunctions(tsComp, ...
    wsComp, w_tsComp, w_ttsComp, epsilon);

% Find where turnover point reaches 1
dsComp = MovingSubstrateFunctions.d(tsComp);
tIdxMaxComp = sum(dsComp <= 1);

% Restrict solutions temporally
tsComp = tsComp(1 : tIdxMaxComp);
wsComp = wsComp(1 : tIdxMaxComp);
w_tsComp = w_tsComp(1 : tIdxMaxComp);
w_ttsComp = w_ttsComp(1 : tIdxMaxComp);
dsComp = dsComp(1 : tIdxMaxComp);

%% Plot pressures in time
% Select timesteps to plot
tsPlot = 0.025 : 0.05 : 0.3
timesteps = floor((tsPlot + IMPACT_TIME) / DELTA_T);

% Variable colors
colorFreq = floor((length(cmap) / 3) / length(timesteps));

% Set up tiled figure
tileFig = tiledlayout(1, 2, 'Tilespacing', 'compact', 'padding', 'compact');

width = 6;
height = 4.5;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

for timestepIdx = 1 : length(timesteps)
    timestep = timesteps(timestepIdx);
    
    tAnalytical = DELTA_T * timestep - IMPACT_TIME; % Analytical time

    % Set line colors
    AnalyticalLineColor = cmap(length(cmap) - (timestepIdx - 1) * colorFreq, :);
    DNSLineColor = cmap((timestepIdx - 1) * colorFreq + 1, :);

    %% Stationary solution for pressure
    nexttile(1);
    hold on;

    % Plot DNS solutions
    plate_output_mat = readmatrix(sprintf("%s/cleaned_data/plate_outputs/output_%d.txt", stat_dir, timestep));
    xs = plate_output_mat(:, 1);
    [sort_xs_stat, sort_idxs] = sort(xs);
    ps = plate_output_mat(sort_idxs, 3);
    stressesStat = plate_output_mat(sort_idxs, 4);

    a = plot(sort_xs_stat, ps, 'color', DNSLineColor, 'linewidth', lineWidth);

    if timestepIdx == 1
        h(1) = a;
    end

    % Plot analytical solution
    [psMoving, ~, ~] ...
            = substratepressure(xsAnalytical, tAnalytical, ...
            StatSubstrateFunctions);

    b = plot(xsAnalytical, psMoving, 'linestyle', ':', ...
        'color', AnalyticalLineColor, 'linewidth', 1.25 * lineWidth);
    if timestepIdx == 1
        h(2) = b;
    end

    xlim([0, 1.25]);
    xticks(0 : 0.25 : 1.5);
    grid on;
    box on;
    xlabel("$r$");
    ylabel("$p(r, 0, t)$");

    title("(a) Stationary plate.", 'Fontsize', fontsize);
    set(gca, 'TitleFontSizeMultiplier', 1);

    %% Moving solution for pressure
    nexttile(2);
    hold on;
    % Plot DNS solutions
    plate_output_mat = readmatrix(sprintf("%s/cleaned_data/plate_outputs/output_%d.txt", moving_dir, timestep));
    xs = plate_output_mat(:, 1);
    [sort_xs_moving, sort_idxs] = sort(xs);
    ps = plate_output_mat(sort_idxs, 3);
    ps(1) = nan;
    stressesMoving = plate_output_mat(sort_idxs, 4);

    plot(sort_xs_moving, ps, 'color', DNSLineColor, 'linewidth', lineWidth);

    % Plot analytical solution
    [psMoving, ~, ~] ...
            = substratepressure(xsAnalytical, tAnalytical, ...
            MovingSubstrateFunctions);
    plot(xsAnalytical, psMoving, 'linestyle', ':', ...
        'color', AnalyticalLineColor, 'linewidth', 1.25 * lineWidth);
    
    xlim([0, 1.25]);
    xticks(0 : 0.25 : 1.5);
    grid on;
    box on;
    xlabel("$r$");
    ylabel("$p(r, -w(t), t)$");

    title("(b) Moving plate.", 'Fontsize', fontsize);
    set(gca, 'TitleFontSizeMultiplier', 1);

    %% Solution for the stresses
%     nexttile(2);
%     hold on;
%     plot(-sort_xs_stat, stressesStat, 'color', blueCol, 'linewidth', lineWidth);
%     plot(sort_xs_moving, stressesMoving, 'color', blueCol, 'linewidth', lineWidth);


end

%% Set figure options
lh = legend([h(1), h(2)], ...
    ["DNS", "Analytical"], ...
    'NumColumns', 2);
lh.Layout.Tile = 'North'; 

set(gcf, 'Renderer', 'Painters');
pause(0.5);

figname = "PlateFigures/PressureComparison";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
% saveas(gcf, sprintf("%s.png", figname));










