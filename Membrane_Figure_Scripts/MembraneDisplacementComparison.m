%% MembraneDisplacementComparison
% 

clear;
close all;

% Add paths
addpath("../Analytical_Scripts/MembraneSolution/FiniteDifference/");
addpath("../Analytical_Scripts/MembraneSolution/FiniteDifference/PressuresFD/");
addpath("../Analytical_Scripts/MembraneSolution/NormalModes/");

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

%% Data directory definitions
% Parent directory where data is stored
parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";

param_dir = append(parent_dir, "RubberSheet");

% DNS directory
dns_dir = "/home/michael/scratch/DPhil_DNS_Data/RubberRuns/ALPHA_1.1_GAMMA_668";

% Stationary directories
stat_FD_dir = "/home/michael/scratch/AnalyticalMembraneTests/Stationary";

%% Parameters
EPSILON = 1;
L = 16;
IMPACT_TIME = 0.125;
T_MAX = 0.35;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 21848;

xs = linspace(0, L, N_MEMBRANE)';

ts_analytical = 0 : DELTA_T : T_MAX;
ts = -IMPACT_TIME : DELTA_T : T_MAX;
impact_timestep = length(ts) - length(ts_analytical) + 1

%% Load normal modes solutions
SolStruct = load(sprintf("%s/NormalModes/SolStruct.mat", param_dir)).SolStruct;
as = SolStruct.as;
a_ts = SolStruct.a_ts;
q_ts = SolStruct.q_ts;
ds = SolStruct.ds;
N_M = SolStruct.N;

%% Plot membrane in time
tPlots = [0, 0.01, 0.1];
timesteps = impact_timestep + tPlots / DELTA_T;


xTicks = 0 : 4 : L;
tiledlayout(1, 2);
width = 6;
height = 2;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

for timestepIdx = 1 : length(timesteps)
    k = timesteps(timestepIdx)
    t = ts(k)

    %% Normal modes solution
    [ws, w_ts, ps] ...
        = MembraneSolutionNM(xs, as(k - impact_timestep + 1, :), ...
        a_ts(k - impact_timestep + 1, :), ...
        q_ts(k - impact_timestep + 1, :), ...
        ds(k - impact_timestep + 1), L, N_M, EPSILON);

    %% Finite difference solution
    ws_FD = load(sprintf("%s/FiniteDifference/composite/w_%d.mat", param_dir, k - impact_timestep)).w_next;
    
    w_ts_FD = load(sprintf("%s/FiniteDifference/composite/w_t_%d.mat", param_dir, k - impact_timestep)).w_t;
    ps_FD = load(sprintf("%s/FiniteDifference/composite/p_%d.mat", param_dir, k - impact_timestep)).p;

    %% DNS solution
    dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", dns_dir, k - 1);
    dns_mat = importdata(dns_filename);
    xs_dns = dns_mat(:, 1);
    [sort_xs, sort_idxs] = sort(xs_dns);
    ws_dns = dns_mat(sort_idxs, 2);
    w_ts_dns = dns_mat(sort_idxs, 3);

    %% Plot displacement
    nexttile(1);
    h(1) = plot(sort_xs, -ws_dns, 'Color', blueCol, ...
        'LineWidth', lineWidth, 'Displayname', 'DNS');
    hold on;
    h(2) = plot(xs, -ws, 'Color', redCol, 'LineStyle', '--', ...
        'LineWidth', lineWidth, 'DisplayName', 'Analytical: NM');
    h(3) = plot(xs(1 : end - 1), -ws_FD, 'Color', redCol, ...
        'LineStyle', ':', 'LineWidth', 1.25 *lineWidth, ...
        'DisplayName', 'Analytical: FD');
    hold off;
    xlim([0, L]);
    xticks(xTicks);
    ylim('padded');
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$-w(x, t)$");

    %% Plot velocity
    nexttile(2);
    plot(sort_xs, -w_ts_dns, 'Color', blueCol, ...
        'LineWidth', lineWidth);
    hold on;
    plot(xs, -w_ts, 'Color', redCol, 'LineStyle', '--', ...
        'LineWidth', lineWidth);
    size(xs)
    size(w_ts_FD)
    plot(xs(1 : end - 1), -w_ts_FD, 'Color', redCol, 'LineStyle', ':', ...
        'LineWidth', 1.25 * lineWidth);
    hold off;
    xlim([0, L]);
    xticks(xTicks);
    ylim('padded');
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$-w_t(x, t)$");

    %% Legend
    % Set figure options
%     lh = legend(h(1 : 3), ...
%         'NumColumns', 3);
%     lh.Layout.Tile = 'North'; 
    
    %% Draw
    drawnow;
    pause(1);

    %% Output
    figname = append("MembraneFigures/Membrane_", num2str(t));
    exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
    
end

%% Plot pressure in time
tPlots = 0.01 : 2e-2 : 0.12;
timesteps = impact_timestep + tPlots / DELTA_T;

xTicks = 0 : 4 : L;
% tiledlayout(2, 1);
figure(2)
width = 6;
height = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Variable colors
colorFreq = floor((length(cmap) / 3) / length(timesteps));

for timestepIdx = 1 : length(timesteps)
    k = timesteps(timestepIdx)
    t = ts(k)

    % Set line colors
    AnalyticalLineColor = cmap(length(cmap) - (timestepIdx - 1) * colorFreq, :);
    DNSLineColor = cmap((timestepIdx - 1) * colorFreq + 1, :);

    % Finite difference solution
    ps_FD = load(sprintf("%s/FiniteDifference/composite/p_%d.mat", param_dir, k - impact_timestep)).p;

    % DNS solution
    dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", dns_dir, k - 1);
    dns_mat = importdata(dns_filename);
    xs_dns = dns_mat(:, 1);
    [sort_xs, sort_idxs] = sort(xs_dns);
    ps_dns = dns_mat(sort_idxs, 4);

    % Plot pressure
    hold on;
    if timestepIdx == 1
        h(1) = plot(sort_xs, ps_dns, 'Color', DNSLineColor, ...
            'LineWidth', lineWidth, 'Displayname', 'DNS');
        h(2) = plot(xs(1 : end - 1), ps_FD, 'Color', AnalyticalLineColor, ...
            'LineStyle', ':', 'LineWidth', 1.25 *lineWidth, ...
            'DisplayName', 'Analytical: FD');
    else
        plot(sort_xs, ps_dns, 'Color', DNSLineColor, ...
            'LineWidth', lineWidth, 'Displayname', 'DNS');
        plot(xs(1 : end - 1), ps_FD, 'Color', AnalyticalLineColor, ...
            'LineStyle', ':', 'LineWidth', 1.25 *lineWidth, ...
            'DisplayName', 'Analytical: FD');
    end
    xlim([0, 0.8]);
    ylim('padded');
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$p(x, -w(x, t), t)$");
end

% Legend
lh = legend(h(1 : 2), 'NumColumns', 3, 'Location', 'NorthEast');

% Output
figname = "MembraneFigures/MembranePressure"
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
    

%% Plot turnovers
ds_FD = load(sprintf("%s/FiniteDifference/composite/ds.mat", param_dir)).ds;
ds_stat = 2 * sqrt(ts_analytical);

% Find dIdx
dIdxNM = sum(ds < 1);
dIdxFD = sum(ds_FD < 1);
dIdxStat = sum(ds_stat < 1);

close(figure(2))
figure(2);
hold on;
plot(ts_analytical(1 : dIdxStat), ds_stat(1 : dIdxStat));
plot(ts_analytical(1 : dIdxNM), ds(1 : dIdxNM));
plot(ts_analytical(1 : dIdxFD), ds_FD(1 : dIdxFD));
