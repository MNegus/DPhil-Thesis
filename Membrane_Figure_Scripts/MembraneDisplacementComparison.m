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

% Stationary DNS directory

stat_dns_dir = "/home/michael/scratch/DPhil_DNS_Data/Stationary_Membrane";

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
N_M = SolStruct.N
d_ts = SolStruct.d_ts;
Js = SolStruct.Js;
Hs = (1 + 4 / pi) * Js;
E_outers = SolStruct.E_outers;
E_jets = SolStruct.E_jets;

%% Plot normal modes stuff
figure(17);
plot(ts_analytical, Js)

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
    pause(0.1);

    %% Output
    figname = append("MembraneFigures/Membrane_", num2str(t));
    exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
    
end

%% Plot pressure in time
tPlots = 0.01 : 2e-2 : 0.12;
timesteps = impact_timestep + tPlots / DELTA_T;

xTicks = 0 : 4 : L;
tiledlayout(2, 1);
% figure(2)
width = 6;
height = 5;
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

    % Normal modes solution
    [~, ~, ps] ...
        = MembraneSolutionNM(xs, as(k - impact_timestep + 1, :), ...
        a_ts(k - impact_timestep + 1, :), ...
        q_ts(k - impact_timestep + 1, :), ...
        ds(k - impact_timestep + 1), L, N_M, EPSILON);

    % Adjusts normal modes past turnover point
    d = ds(k - impact_timestep + 1);
    psEndIdx = sum(xs < d);
    ps(xs >= d) = nan;

    % Finite difference solution
    ps_FD = load(sprintf("%s/FiniteDifference/composite/p_%d.mat", param_dir, k - impact_timestep)).p;

    % DNS solution
    dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", dns_dir, k - 1);
    dns_mat = importdata(dns_filename);
    xs_dns = dns_mat(:, 1);
    [sort_xs, sort_idxs] = sort(xs_dns);
    ps_dns = dns_mat(sort_idxs, 4);

    % Plot pressure
    if timestepIdx == 1
        nexttile(1);
        hold on;
        h(1) = plot(sort_xs, ps_dns, 'Color', DNSLineColor, ...
            'LineWidth', lineWidth, 'Displayname', 'DNS');
        h(2) = plot(xs(1 : end - 1), ps_FD, 'Color', AnalyticalLineColor, ...
            'LineStyle', ':', 'LineWidth', 1.25 *lineWidth, ...
            'DisplayName', 'Analytical: FD');

        nexttile(2);
        hold on;
        h(3) = plot(xs, ps, 'Color', AnalyticalLineColor, 'LineStyle', '--', ...
            'LineWidth', lineWidth, 'Displayname', 'Analytical: NM');
        line([d d], [0 ps(psEndIdx)], 'Color', AnalyticalLineColor, 'LineStyle', '--');
    else
        nexttile(1);
        hold on;
        plot(sort_xs, ps_dns, 'Color', DNSLineColor, ...
            'LineWidth', lineWidth, 'Displayname', 'DNS');
        plot(xs(1 : end - 1), ps_FD, 'Color', AnalyticalLineColor, ...
            'LineStyle', ':', 'LineWidth', 1.25 *lineWidth, ...
            'DisplayName', 'Analytical: FD');

        nexttile(2);
        hold on;
        plot(xs, ps, 'Color', AnalyticalLineColor, 'LineStyle', '--', ...
            'LineWidth', lineWidth, 'Displayname', 'Analytical: NM');
        
        line([d d], [0 ps(psEndIdx)], 'Color', AnalyticalLineColor, 'LineStyle', '--');
    end

    nexttile(1);
    xlim([0, 0.8]);
    ylim([-2, 50]);
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$p(x, -w(x, t), t)$");

    nexttile(2);
    xlim([0, 0.8]);
    ylim([-2, 20]);
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$p(x, -w(x, t), t)$");
end

% Legend
lh = legend([h(1), h(3), h(2)], 'NumColumns', 3);
lh.Layout.Tile = 'North'; 

% Output
figname = "MembraneFigures/MembranePressure"
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
    

%% Turnover and jet thickness comparison

% Time axes settings
xLims = [0, 0.3];

% Load stationary DNS
turnover_mat = importdata(sprintf("%s/turnover_points.csv", stat_dns_dir));
ts_stat_DNS = turnover_mat(:, 1) - IMPACT_TIME;
ds_stat_DNS = turnover_mat(:, 2);
Hs_stat_DNS = turnover_mat(:, 3);

% Load stationary analytical
dsStat = 2 * sqrt(ts_analytical);
d_tStat = 1 ./ sqrt(ts_analytical);
JsStat = pi * ts_analytical.^(3/2) / 4;
HsStat = (1 + 4 / pi) * JsStat;

% Load moving DNS
turnover_mat = importdata(sprintf("%s/turnover_points.csv", dns_dir));
ts_DNS = turnover_mat(:, 1) - IMPACT_TIME;
ds_DNS = turnover_mat(:, 2);
Hs_DNS = turnover_mat(:, 3);

% Load moving FD
ds_FD = load(sprintf("%s/FiniteDifference/composite/ds.mat", param_dir)).ds;
d_ts_FD = load(sprintf("%s/FiniteDifference/composite/d_ts.mat", param_dir)).d_ts;
Js_FD = load(sprintf("%s/FiniteDifference/composite/Js.mat", param_dir)).Js;
Hs_FD = (1 + 4 / pi) * Js_FD;

% Find index of max turnover times
dIdxStat = sum(dsStat < 1);
dIdxNM = sum(ds < 1);
dIdxFD = sum(ds_FD < 1);

%% Turnover point plot

tiledlayout(1, 2);
width = 6;
height = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Stationary membrane solution
nexttile;
hold on;
xline(ts_analytical(dIdxStat), '--', 'color', 0.75 * [1 1 1]);
xline(ts_analytical(dIdxStat), ':', 'color', 0.75 * [1 1 1]);
h(1) = plot(ts_stat_DNS, ds_stat_DNS, 'color', blueCol, 'linewidth', lineWidth);
h(2) = plot(ts_analytical(1 : dIdxStat), dsStat(1 : dIdxStat), 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
h(3) = plot(ts_analytical(1 : dIdxStat), dsStat(1 : dIdxStat), 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);
grid on;
box on;
xlabel("$t$");
ylabel("$d(t)$");
xlim(xLims);
xticks([-0.1 : 0.1 : 0.3]);
ylim([0, 1.2]);
title("(a) Stationary membrane.", 'Fontsize', fontsize);
set(gca, 'TitleFontSizeMultiplier', 1);

% Moving membrane solution
nexttile;
hold on;
xline(ts_analytical(dIdxNM), '--', 'color', 0.75 * [1 1 1]);
xline(ts_analytical(dIdxFD), ':', 'color', 0.75 * [1 1 1]);
h(1) = plot(ts_DNS, ds_DNS, 'color', blueCol, 'linewidth', lineWidth);
h(2) = plot(ts_analytical(1 : dIdxNM), ds(1 : dIdxNM), 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
h(3) = plot(ts_analytical(1 : dIdxFD), ds_FD(1 : dIdxFD), 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);
grid on;
box on;
xlabel("$t$");
ylabel("$d(t)$");
xlim(xLims);
xticks([-0.1 : 0.1 : 0.3]);
ylim([0, 1.2]);
title("(b) Moving membrane.", 'Fontsize', fontsize);
set(gca, 'TitleFontSizeMultiplier', 1);

% Set figure options
lh = legend(h(1 : 3), ...
    ["DNS", "Analytical: NM", "Analytical: FD"], ...
    'NumColumns', 3);
lh.Layout.Tile = 'North'; 

figname = "MembraneFigures/MembraneTurnoverPointComparison";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Jet height
tiledlayout(1, 2);
width = 6;
height = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Stationary membrane solution
nexttile;
hold on;
xline(ts_analytical(dIdxStat), '--', 'color', 0.75 * [1 1 1]);
xline(ts_analytical(dIdxStat), ':', 'color', 0.75 * [1 1 1]);
h(1) = plot(ts_stat_DNS, Hs_stat_DNS, 'color', blueCol, 'linewidth', lineWidth);
h(2) = plot(ts_analytical(1 : dIdxStat), HsStat(1 : dIdxStat), 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
h(3) = plot(ts_analytical(1 : dIdxStat), HsStat(1 : dIdxStat), 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);
grid on;
box on;
xlabel("$t$");
xlim(xLims);
xticks([-0.1 : 0.1 : 0.3]);
ylabel("$H(t)$");
title("(a) Stationary membrane.", 'Fontsize', fontsize);
set(gca, 'TitleFontSizeMultiplier', 1);


% Moving membrane solution
nexttile;
hold on;
xline(ts_analytical(dIdxNM), '--', 'color', 0.75 * [1 1 1]);
xline(ts_analytical(dIdxFD), ':', 'color', 0.75 * [1 1 1]);
plot(ts_DNS, Hs_DNS, 'color', blueCol, 'linewidth', lineWidth);
plot(ts_analytical(1 : dIdxNM), Hs(1 : dIdxNM), 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
plot(ts_analytical(1 : dIdxFD), Hs_FD(1 : dIdxFD), 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);
grid on;
box on;
xlabel("$t$");
xlim(xLims);
xticks([-0.1 : 0.1 : 0.3]);
ylabel("$H(t)$");
title("(b) Moving membrane.", 'Fontsize', fontsize);
set(gca, 'TitleFontSizeMultiplier', 1);

% Set figure options
lh = legend(h(1 : 3), ...
    ["DNS", "Analytical: NM", "Analytical: FD"], ...
    'NumColumns', 3);
lh.Layout.Tile = 'North'; 

set(gcf, 'Renderer', 'Painters');
pause(0.5);

figname = "MembraneFigures/MembraneTurnoverHeightComparison";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

(HsStat(dIdxStat) - Hs(dIdxStat)) / HsStat(dIdxStat)

(Hs_stat_DNS(375) - Hs_DNS(375)) / Hs_stat_DNS(375)

%% Maximum pressure
% DNS times
timesteps = 0 : 7000;
ts_dns = timesteps * DELTA_T - IMPACT_TIME;

% Load stationary DNS solutions
StatMaxStruct = load(sprintf("%s/MaxStruct.mat", stat_dns_dir)).MaxStruct;
psMaxs_stat_DNS = StatMaxStruct.pMaxs;

% Load moving DNS solutions
MovingMaxStruct = load(sprintf("%s/MaxStruct.mat", dns_dir)).MaxStruct;
psMaxs_DNS = MovingMaxStruct.pMaxs;

% Load FD solutions
pMaxStat = d_tStat.^2 / 2; 
pMaxs_FD = d_ts_FD.^2 / 2;

% Load NM solutions
pMaxs_NM = d_ts.^2 / 2;

% Tiled plot
% tiledlayout(1, 2);
close(figure(46));
figure(46);
hold on;
width = 6;
height = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

sz = 15;
hold on;
scatter(ts_dns, psMaxs_DNS, sz, cmap(100, :));
plot(ts_analytical, pMaxs_NM, 'linestyle', '--', ...
    'color', redCol, 'linewidth', lineWidth);
plot(ts_analytical(2 : end - 1), pMaxs_FD(2 : end - 1), 'linestyle', ':', ...
    'color', redCol, 'linewidth', 1.25 * lineWidth);
ylim([0, 10]);
xlim([-0.05, max(ts_dns)])
grid on;
xlabel("$t$");
ylabel("max$(p(x, t))$");

lh = legend("DNS", "Analytical: NM", "Analytical: FD", 'Location', 'northeast');

figname = "MembraneFigures/MaxPressureComparison";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Mass loss plot
stat_output_mat = importdata(sprintf("%s/output.txt", stat_dns_dir));
stat_dns_ts = stat_output_mat(:, 1) - IMPACT_TIME;
stat_dns_area = stat_output_mat(:, 2) / (2 * pi);
stat_mass_ratio = (stat_dns_area - pi / 2) / (pi / 2);

output_mat = importdata(sprintf("%s/output.txt", dns_dir));
dns_ts = output_mat(:, 1) - IMPACT_TIME;
dns_area = output_mat(:, 2) / (2 * pi);
mass_ratio = (dns_area - pi / 2) / (pi / 2);

close(figure(56));
figure(56);
hold on;
width = 5;
height = 2;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);
scatter(stat_dns_ts, stat_mass_ratio, sz, blueCol);
scatter(dns_ts, mass_ratio, sz, redCol);
ylim("padded");
xlim("padded");
grid on;
box on;
xlabel("$t$");
ylabel("$R(t)$");
legend("Stationary membrane", "Moving membrane", 'Location', 'southwest', 'Numcolumns', 1);

figname = "MembraneFigures/MembraneMassLoss";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Plot energy
close(figure(65));
figure(65);
hold on;
plot(ts_analytical, E_outers);
plot(ts_analytical, E_jets);
plot(ts_analytical, pi * ts_analytical);


%% Turnover point (STATIONARY COMPARISON)
% 
% % Time axes settings
% xLims = [0, 0.32];
% 
% % Load stationary DNS
% turnover_mat = importdata(sprintf("%s/turnover_points.csv", stat_dns_dir))
% ts_DNS_stat = turnover_mat(:, 1) - IMPACT_TIME;
% ds_DNS_stat = turnover_mat(:, 2);
% 
% % Load stationary analytical
% dsStat = 2 * sqrt(ts_analytical);
% 
% % Load moving DNS
% turnover_mat = importdata(sprintf("%s/turnover_points.csv", dns_dir));
% ts_DNS = turnover_mat(:, 1) - IMPACT_TIME;
% ds_DNS = turnover_mat(:, 2);
% 
% % Load moving FD
% ds_FD = load(sprintf("%s/FiniteDifference/composite/ds.mat", param_dir)).ds;
% d_ts_FD = load(sprintf("%s/FiniteDifference/composite/d_ts.mat", param_dir)).d_ts;
% 
% % Find index of max turnover times
% dIdxStat = sum(dsStat < 1);
% dIdxNM = sum(ds < 1);
% dIdxFD = sum(ds_FD < 1);
% 
% tiledlayout(1, 2);
% width = 6;
% height = 3;
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);
% 
% % Stationary solutions
% nexttile;
% hold on;
% xline(ts_analytical(dIdxStat), '--', 'color', 0.75 * [1 1 1]);
% xline(ts_analytical(dIdxStat), ':', 'color', 0.75 * [1 1 1]);
% h(1) = plot(ts_DNS_stat, ds_DNS_stat, 'color', blueCol, 'linewidth', lineWidth);
% h(2) = plot(ts_analytical(1 : dIdxStat), dsStat(1 : dIdxStat), 'linestyle', '--', ...
%     'color', redCol, 'linewidth', lineWidth);
% h(3) = plot(ts_analytical(1 : dIdxStat), dsStat(1 : dIdxStat), 'linestyle', ':', ...
%     'color', redCol, 'linewidth', 1.25 * lineWidth);
% 
% grid on;
% box on;
% 
% xlabel("$t$");
% xlim(xLims);
% xticks([-0.1 : 0.1 : 0.4]);
% 
% ylabel("$d(t)$");
% ylim([0, 1.2]);
% 
% title("(a) Stationary membrane.", 'Fontsize', fontsize);
% set(gca, 'TitleFontSizeMultiplier', 1);
% 
% % Moving solutions
% nexttile;
% hold on;
% xline(ts_analytical(dIdxNM), '--', 'color', 0.75 * [1 1 1]);
% xline(ts_analytical(dIdxFD), ':', 'color', 0.75 * [1 1 1]);
% plot(ts_DNS, ds_DNS, 'color', blueCol, 'linewidth', lineWidth);
% plot(ts_analytical(1 : dIdxNM), ds(1 : dIdxNM), 'linestyle', '--', ...
%     'color', redCol, 'linewidth', lineWidth);
% plot(ts_analytical(1 : dIdxFD), ds_FD(1 : dIdxFD), 'linestyle', ':', ...
%     'color', redCol, 'linewidth', 1.25 * lineWidth);
% 
% grid on;
% box on;
% 
% xlabel("$t$");
% xlim(xLims);
% xticks([-0.1 : 0.1 : 0.4]);
% 
% ylabel("$d(t)$");
% ylim([0, 1.2]);
% 
% title("(b) Moving membrane.", 'Fontsize', fontsize);
% set(gca, 'TitleFontSizeMultiplier', 1);
% 
% 
% % Set figure options
% % lh = legend(h(1 : 3), ...
% %     ["DNS", "Analytical (leading-order)", "Analytical (composite)"], ...
% %     'NumColumns', 3);
% % lh.Layout.Tile = 'North'; 
% 
% set(gcf, 'Renderer', 'Painters');
% pause(0.5);
% 
% figname = "MembraneFigures/MembraneTurnoverComparison";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
