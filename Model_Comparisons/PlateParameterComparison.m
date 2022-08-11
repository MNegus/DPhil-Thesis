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

%% Parameter varying values
% ALPHA varying parameters
ALPHAStruct.ALPHA_strs =  ["2.0", "5.0", "10.0", "20.0", "100.0"];
ALPHAStruct.BETA_strs = ["0.0", "0.0", "0.0", "0.0", "0.0"];
ALPHAStruct.GAMMA_strs = ["0.0", "0.0", "0.0", "0.0", "0.0"];
ALPHAStruct.noElements = length(ALPHAStruct.ALPHA_strs);
ALPHAStruct.type = "ALPHA";

% BETA varying parameters
BETAStruct.BETA_strs = ["0.0", "7.07", "28.28", "141.42"];
BETAStruct.ALPHA_strs = ["2.0", "2.0", "2.0", "2.0"];
BETAStruct.GAMMA_strs = ["100.0", "100.0", "100.0", "100.0"];
BETAStruct.noElements = length(BETAStruct.BETA_strs);
BETAStruct.type = "BETA";

% GAMMA varying parameters
GAMMAStruct.GAMMA_strs = ["10.0", "20.0", "100.0", "500.0", "1000.0"];
GAMMAStruct.ALPHA_strs = ["2.0", "2.0", "2.0", "2.0", "2.0"];
GAMMAStruct.BETA_strs = ["0.0", "0.0", "0.0", "0.0", "0.0"];
GAMMAStruct.noElements = length(GAMMAStruct.GAMMA_strs);
GAMMAStruct.type = "GAMMA";

% Combine structs
ParameterStructs = [ALPHAStruct, BETAStruct, GAMMAStruct];


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

%% Loop over varying types
for ParameterStruct = ParameterStructs
    % Load parameters
    ALPHA_strs = ParameterStruct.ALPHA_strs;
    BETA_strs = ParameterStruct.BETA_strs;
    GAMMA_strs = ParameterStruct.GAMMA_strs;
    noElements = ParameterStruct.noElements;
    type = ParameterStruct.type;

    % Set up colors
    extent = 3;
    colorFreq = floor((length(cmap) / extent) / noElements);
    
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
    
    % Plot parameter varying DNS solutions
    for ParamIdx = 1 : noElements

        % Set type
        if type == "ALPHA"
            paramStr = ALPHA_strs(ParamIdx);
        elseif type == "BETA"
            paramStr = BETA_strs(ParamIdx);
        else
            paramStr = GAMMA_strs(ParamIdx);
        end
    
        % Set line colors
        DNSLineColor = cmap(DNSLineColorIdxs(ParamIdx), :);
    
        % Load DNS solutions
        param_dir = append(dns_master_dir, ...
            "/Plate_Parameter_Runs/", type, "_varying/", type, "_", paramStr);
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
    for ParamIdx = 1 : noElements
        % Set type
        if type == "ALPHA"
            paramStr = ALPHA_strs(ParamIdx);
        elseif type == "BETA"
            paramStr = BETA_strs(ParamIdx);
        else
            paramStr = GAMMA_strs(ParamIdx);
        end
    
        % Set line colors
        AnalyticalLineColor = cmap(AnalyticalLineColorIdxs(ParamIdx), :);
    
        % Load analytical solutions struct
        fileName = append("AnalyticalSolutions/", type, "_varying/" , type, "_", paramStr, ".mat")
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
    mkdir(append("PlateFigures/", type, "_varying"));
    figname = append("PlateFigures/", type, "_varying/Plate_Force_", type);
    exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
end

% %% ALPHA varying
% BETA = 0;
% GAMMA = 0;
% ALPHA_strs = ["2.0", "5.0", "10.0", "20.0", "100.0"];
% 
% % Set up colors
% extent = 3;
% colorFreq = floor((length(cmap) / extent) / length(ALPHA_strs));
% 
% % Line color indices
% DNSLineColorIdxs = floor(length(cmap) / extent) : -colorFreq : 1;
% AnalyticalLineColorIdxs = (length(cmap) - floor(length(cmap) / extent)) : colorFreq : length(cmap);
% 
% %% Substrate and force comparison
% tiledlayout(1, 2);
% width = 6;
% height = 4;
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);
% 
% % Plot stationary DNS solution
% nexttile(2);
% hold on;
% h(1) = plot(ts, FsStat, 'Color', 'black', 'LineWidth', lineWidth);
% 
% % Plot alpha varying DNS solutions
% for ALPHA_idx = 1 : length(ALPHA_strs)
%     ALPHA_str = ALPHA_strs(ALPHA_idx);
%     
%     % Load numerical value for ALPHA
%     ALPHA = str2double(ALPHA_str);
% 
%     % Set line colors
%     DNSLineColor = cmap(DNSLineColorIdxs(ALPHA_idx), :);
% 
%     % Load DNS solutions
%     param_dir = append(dns_master_dir, ...
%         "/Plate_Parameter_Runs/ALPHA_varying/ALPHA_", ALPHA_str);
%     output_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", param_dir));
%     ts = output_mat(:, 1) - IMPACT_TIME;
%     FsDNS = output_mat(:, 3);
%     wsDNS = output_mat(:, 6);
%     w_tsDNS = output_mat(:, 7);
%     w_ttsDNS = output_mat(:, 8);
% 
%     % Plot substrate deformation and force
%     nexttile(1);
%     hold on;
%     plot(ts, -wsDNS, 'Color', DNSLineColor, 'LineWidth', lineWidth);
% 
%     % Force
%     nexttile(2);
%     hold on;
%     plot(ts, FsDNS, 'Color', DNSLineColor, 'LineWidth', lineWidth);
% 
% 
% end
% 
% % Plot alpha varying analytical solutions
% for ALPHA_idx = 1 : length(ALPHA_strs)
%     ALPHA_str = ALPHA_strs(ALPHA_idx)
% 
%     % Set line colors
%     AnalyticalLineColor = cmap(AnalyticalLineColorIdxs(ALPHA_idx), :);
% 
%     % Load analytical solutions struct
%     fileName = append("AnalyticalSolutions/ALPHA_varying/ALPHA_", ALPHA_str, ".mat");
%     SolStruct = load(fileName).SolStruct;
% 
%     % Plot substrate deformation and force
%     nexttile(1);
%     hold on;
%     plot(SolStruct.ts, -SolStruct.ws, 'Color', AnalyticalLineColor, ...
%         'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);
% 
%     nexttile(2);
%     hold on;
%     plot(SolStruct.ts, SolStruct.FsComp, 'Color', AnalyticalLineColor, ...
%         'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);
% 
% end
% 
% % Plot stationary analytical solution
% h(2) = plot(StatStruct.ts, StatStruct.FsComp, 'Color', 'black', ...
%         'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);
% 
% % Set figure options
% nexttile(1);
% grid on;
% box on;
% ylim("padded")
% xlim([-0.2, 0.7]);
% xlabel("$t$");
% ylabel("$-w(t)$");
% 
% nexttile(2);
% grid on;
% box on;
% ylim("padded");
% xlim([-0.2, 0.7]);
% xlabel("$t$");
% ylabel("$F(t)$");
% 
% lh = legend(h(1 : 2), ["DNS", "Analytical"], 'NumColumns', 2);
% lh.Layout.Tile = 'North'; 
% 
% % Export figure
% mkdir("PlateFigures/ALPHA_varying");
% figname = "PlateFigures/ALPHA_varying/Plate_Force_ALPHA";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
% 
% %% Turnover compare
% tMax = 0.45; 
% dsStatTurnover(tsStatTurnover > tMax) = nan;
% JsStatTurnover(tsStatTurnover > tMax) = nan;
% 
% dsStatTurnover(tsStatTurnover < 0) = nan;
% JsStatTurnover(tsStatTurnover < 0) = nan;
% 
% tiledlayout(1, 2);
% width = 6;
% height = 4;
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);
% 
% % Plot alpha varying DNS solutions
% for ALPHA_idx = 1 : length(ALPHA_strs)
%     ALPHA_str = ALPHA_strs(ALPHA_idx);
%     
%     % Load numerical value for ALPHA
%     ALPHA = str2double(ALPHA_str);
% 
%     % Set line colors
%     DNSLineColor = cmap(DNSLineColorIdxs(ALPHA_idx), :);
% 
%     %% Load DNS solutions
%     param_dir = append(dns_master_dir, ...
%         "/Plate_Parameter_Runs/ALPHA_varying/ALPHA_", ALPHA_str);
%     turnover_mat = readmatrix(sprintf("%s/turnover_points.csv", param_dir));
%     
%     % Load solutions
%     ts = turnover_mat(:, 1) - IMPACT_TIME;
%     ds = turnover_mat(:, 2);
%     Hs = turnover_mat(:, 3);
%     Hs(ts < 0) = nan;
%     ds(ts < 0) = nan;
% 
%     % Restrict times
%     ds(ts > tMax) = nan;
%     Hs(ts > tMax) = nan;
% 
%     % Plot turnover point radial position
%     nexttile(1);
%     hold on;
%     plot(ts, ds, 'Color', DNSLineColor, 'LineWidth', lineWidth);
% 
%     % Plot turnover point vertical position
%     nexttile(2);
%     hold on;
%     plot(ts, Hs, 'Color', DNSLineColor, 'LineWidth', lineWidth);
% 
% 
% end
% 
% % Plot stationary DNS solution
% nexttile(1);
% hold on;
% h(1) = plot(tsStatTurnover, dsStatTurnover, 'Color', 'black', 'LineWidth', lineWidth);
% 
% nexttile(2);
% hold on;
% plot(tsStatTurnover, JsStatTurnover, 'Color', 'black', 'LineWidth', lineWidth);
% 
% % % Plot alpha varying analytical solutions
% % for ALPHA_idx = 1 : length(ALPHA_strs)
% %     ALPHA_str = ALPHA_strs(ALPHA_idx)
% % 
% %     % Set line colors
% %     AnalyticalLineColor = cmap(AnalyticalLineColorIdxs(ALPHA_idx), :);
% % 
% %     % Load analytical solutions struct
% %     fileName = append("AnalyticalSolutions/ALPHA_varying/ALPHA_", ALPHA_str, ".mat");
% %     SolStruct = load(fileName).SolStruct;
% % 
% %     % Plot substrate deformation and force
% %     nexttile(1);
% %     hold on;
% %     plot(SolStruct.ts, SolStruct.ds, 'Color', AnalyticalLineColor, ...
% %         'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);
% % 
% %     nexttile(2);
% %     hold on;
% %     plot(SolStruct.ts, (1 + 4 / pi) * SolStruct.Js, 'Color', AnalyticalLineColor, ...
% %         'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);
% % 
% % end
% % 
% % % Plot stationary analytical solution
% % nexttile(1);
% % hold on;
% % h(2) = plot(StatStruct.ts, StatStruct.ds, 'Color', 'black', ...
% %         'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);
% % 
% % nexttile(2);
% % hold on;
% % plot(StatStruct.ts, (1 + 4 / pi) * StatStruct.Js, 'Color', 'black', ...
% %         'LineStyle', ':', 'LineWidth', 1.25 * lineWidth);
% 
% % Set figure options
% nexttile(1);
% grid on;
% box on;
% ylim("padded")
% xlim([-0.1, 0.5]);
% xlabel("$t$");
% ylabel("$d(t)$");
% 
% nexttile(2);
% grid on;
% box on;
% ylim("padded");
% xlim([-0.1, 0.5]);
% xlabel("$t$");
% ylabel("$H(t)$");
% 
% % lh = legend(h(1), "DNS", 'NumColumns', 2);
% % lh.Layout.Tile = 'North'; 
% 
% % Export figure
% mkdir("PlateFigures/ALPHA_varying");
% figname = "PlateFigures/ALPHA_varying/Turnover_ALPHA";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
% 
