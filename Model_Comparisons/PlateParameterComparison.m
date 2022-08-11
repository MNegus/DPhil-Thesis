%% PlateAndForceComparison.m
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

%% Load stationary plate DNS solutions
stat_dir = sprintf("%s/Stationary_Plate/axi", dns_master_dir);
stat_output_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", stat_dir));
ts = stat_output_mat(:, 1) - IMPACT_TIME;
FsStat = stat_output_mat(:, 3);

% Turnover points
stat_turnover_mat = readmatrix(sprintf("%s/turnover_points.csv", stat_dir));
tsStatTurnover = stat_turnover_mat(:, 1) - IMPACT_TIME;
dsStatTurnover = stat_turnover_mat(:, 2);
JsStatTurnover = stat_turnover_mat(:, 3);
JsStatTurnover(tsStatTurnover < 0) = 0;


%% Load stationary analytical solutions
StatStruct = load("AnalyticalSolutions/StationarySol.mat").SolStruct;

%% ALPHA varying
BETA = 0;
GAMMA = 0;
ALPHA_strs = ["2.0", "5.0", "10.0", "20.0", "100.0"];

% Set up colors
extent = 2.5;
colorFreq = floor((length(cmap) / extent) / length(ALPHA_strs));

% Line color indices
DNSLineColorIdxs = floor(length(cmap) / extent) : -colorFreq : 1;
AnalyticalLineColorIdxs = (length(cmap) - floor(length(cmap) / extent)) : colorFreq : length(cmap);

%% Substrate and force comparison
tiledlayout(1, 2);
width = 6;
height = 4;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Plot stationary DNS solution
nexttile(2);
hold on;
h(1) = plot(ts, FsStat, 'Color', 'black', 'LineWidth', lineWidth);

% Plot alpha varying DNS solutions
for ALPHA_idx = 1 : length(ALPHA_strs)
    ALPHA_str = ALPHA_strs(ALPHA_idx);
    
    % Load numerical value for ALPHA
    ALPHA = str2double(ALPHA_str);

    % Set line colors
%     DNSLineColor = cmap((ALPHA_idx - 1) * colorFreq + 1, :);
    DNSLineColor = cmap(DNSLineColorIdxs(ALPHA_idx), :);

    %% Load DNS solutions
    param_dir = append(dns_master_dir, ...
        "/Plate_Parameter_Runs/ALPHA_varying/ALPHA_", ALPHA_str);
    output_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", param_dir));
    ts = output_mat(:, 1) - IMPACT_TIME;
    FsDNS = output_mat(:, 3);
    wsDNS = output_mat(:, 6);
    w_tsDNS = output_mat(:, 7);
    w_ttsDNS = output_mat(:, 8);

    % Plot substrate deformation and force
    nexttile(1);
    hold on;
    plot(ts, -wsDNS, 'Color', DNSLineColor, 'LineWidth', lineWidth);

    % Force
    nexttile(2);
    hold on;
    plot(ts, FsDNS, 'Color', DNSLineColor, 'LineWidth', lineWidth);


end

% Plot alpha varying analytical solutions
for ALPHA_idx = 1 : length(ALPHA_strs)
    ALPHA_str = ALPHA_strs(ALPHA_idx)

    % Set line colors
%     AnalyticalLineColor = cmap(length(cmap) - (ALPHA_idx - 1) * colorFreq, :);
    AnalyticalLineColor = cmap(AnalyticalLineColorIdxs(ALPHA_idx), :);

    % Load analytical solutions struct
    fileName = append("AnalyticalSolutions/ALPHA_varying/ALPHA_", ALPHA_str, ".mat");
    SolStruct = load(fileName).SolStruct;

    % Plot substrate deformation and force
    nexttile(1);
    hold on;
    plot(SolStruct.ts, -SolStruct.ws, 'Color', AnalyticalLineColor, ...
        'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);

    nexttile(2);
    hold on;
    plot(SolStruct.ts, SolStruct.FsComp, 'Color', AnalyticalLineColor, ...
        'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);

end

% Plot stationary analytical solution
h(2) = plot(StatStruct.ts, StatStruct.FsComp, 'Color', 'black', ...
        'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);

% Set figure options
nexttile(1);
grid on;
box on;
ylim("padded")
xlim([-0.2, 0.7]);
xlabel("$t$");
ylabel("$-w(t)$");

nexttile(2);
grid on;
box on;
ylim("padded");
xlim([-0.2, 0.7]);
xlabel("$t$");
ylabel("$F(t)$");

lh = legend(h(1 : 2), ["DNS", "Analytical"], 'NumColumns', 2);
lh.Layout.Tile = 'North'; 

% Export figure
mkdir("PlateFigures/ALPHA_varying");
figname = "PlateFigures/ALPHA_varying/Plate_Force_ALPHA";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Turnover compare


