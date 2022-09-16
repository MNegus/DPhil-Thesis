%% DNSMembraneCompare.m
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

%% Data directories

% Stationary DNS directory
stat_dns_dir = "/home/michael/scratch/DPhil_DNS_Data/Stationary_Membrane";

% Parent directory where DNS data is stored
parent_dir = "/home/michael/scratch/DPhil_DNS_Data/MembraneParameterRuns/";

% Parameter varying directories
delta_dir = append(parent_dir, "DELTA_varying/");
gamma_dir = append(parent_dir, "GAMMA_varying/");
beta_dir = append(parent_dir, "BETA_varying/");
new_beta_dir = append(parent_dir, "NEW_BETA_varying/");
new_gamma_dir = append(parent_dir, "NEW_GAMMA_varying/");

%% Parameters
EPSILON = 1;
L = 16;
IMPACT_TIME = 0.125;
T_MAX = 0.7 - IMPACT_TIME;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 21848;

xs = linspace(0, L, N_MEMBRANE)';

ts_analytical = 0 : DELTA_T : T_MAX;
ts = -IMPACT_TIME : DELTA_T : T_MAX;
timesteps = 0 : length(ts) - 1;
impact_timestep = length(ts) - length(ts_analytical) + 1;

xTicks = 0 : 4 : L;

%% Delta varying
ALPHA_strs = ["0.1375", "0.3889", "1.1", "3.111", "8.8"];
GAMMA_strs = ["1.305", "29.52", "668.0", "15110.0", "342000.0"];
BETA_strs = ["0", "0", "0", "0", "0"];
DELTA_strs = ["0.125", "0.3536", "1.000", "2.828", "8.000"];
noElements = length(ALPHA_strs);
vary_dir = delta_dir;

% Line colours
DNSLineColorIdxs = floor(linspace(1, length(cmap), noElements));

% Setup figure
tiledlayout(1, 2);
width = 6;
height = 2;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Plotting times
tPlots = [0.1, 0.2, 0.575];
timesteps = cast(impact_timestep + tPlots / DELTA_T, "int32")

% Plot membrane in time
for k = timesteps
    t = ts(k);

    % Loop over cases
    for paramIdx = noElements : -1 : 1
        % Set param dir
        param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", vary_dir, ...
            ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));

        % Set line color
        lineColor = cmap(DNSLineColorIdxs(paramIdx), :);

        % Set display name
        displayName = append("$\delta = $", DELTA_strs(paramIdx));

        % Load dns files
        dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1);
        dns_mat = importdata(dns_filename);
        xs_dns = dns_mat(:, 1);
        [sort_xs, sort_idxs] = sort(xs_dns);
        ws_dns = dns_mat(sort_idxs, 2);
        w_ts_dns = dns_mat(sort_idxs, 3);
        ps_dns = dns_mat(sort_idxs, 4);
        
        % Plot membrane position
        nexttile(1);
        h(paramIdx) = plot(sort_xs, -ws_dns, 'LineWidth', lineWidth, ...
            'Color', lineColor, 'Displayname', displayName);
        hold on;

        % Plot membrane derivative
        nexttile(2);
        plot(sort_xs, -w_ts_dns, 'LineWidth', lineWidth, 'Color', lineColor);
        hold on;
    end

    % Set figure options
    nexttile(1);
    hold off;
    xlim([0, 16]);
    xticks(xTicks);
    ylim([-0.45, 0.1])
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$-w(x, t)$");


    nexttile(2);
    hold off;
    xlim([0, 16]);
    xticks(xTicks);
    ylim([-0.85, 0.375])
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$-w_t(x, t)$");

    % Set legend options
%     lh = legend(h(1 : noElements), ...
%         'NumColumns', 3);
%     lh.Layout.Tile = 'North'; 

    % Pause
    pause(0.01);

    % Output
    figname = append("MembraneFigures/MembraneParamVary/MembraneDELTAVary_", num2str(t));
    exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
end

% Plot pressure on substrate
k = timesteps(1)

close(figure(1));
figure(1);
width = 6;
height = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Loop over cases
for paramIdx = 1 : noElements
    % Set param dir
    param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", delta_dir, ...
        ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));

    % Set line color
    lineColor = cmap(DNSLineColorIdxs(paramIdx), :);

    % Set display name
    displayName = append("$\delta = $", DELTA_strs(paramIdx));

    % Load dns files
    dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1);
    dns_mat = importdata(dns_filename);
    xs_dns = dns_mat(:, 1);
    [sort_xs, sort_idxs] = sort(xs_dns);
    ps_dns = dns_mat(sort_idxs, 4);
    
    % Plot membrane position
    h(paramIdx) = plot(sort_xs, ps_dns, 'LineWidth', lineWidth, ...
        'Color', lineColor, 'Displayname', displayName);
    hold on;

end

% Plot stationary pressure
 dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", stat_dns_dir, k - 1);
dns_mat = importdata(dns_filename);
xs_dns = dns_mat(:, 1);
[sort_xs, sort_idxs] = sort(xs_dns);
ps_dns = dns_mat(sort_idxs, 4);
h(noElements + 1) = plot(sort_xs, ps_dns, 'LineWidth', 1.25 * lineWidth, ...
        'Color', 'black', 'Linestyle', ':','Displayname', 'Stationary');


% Set figure options
hold off;
xlim([0, 1]);
ylim('padded')
grid on;
box on;
xlabel("$x$");
ylabel("$p(x, t)$");


% Set legend options
lh = legend(h(1 : noElements + 1), ...
    'NumColumns', 3, 'Location', 'Northoutside');

% Pause
pause(0.01);

% Output
figname = "MembraneFigures/MembraneParamVary/PressureDELTAVary";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);


%% Gamma varying
% GAMMA_strs = ["1.883", "18.83", "188.3", "1883.0", "18830.0"];
% ALPHA_strs = ["1.1", "1.1", "1.1", "1.1", "1.1"];
% BETA_strs = ["0", "0", "0", "0", "0"];
% noElements = length(GAMMA_strs);
% vary_dir = gamma_dir;
% 
% % Line colours
% DNSLineColorIdxs = floor(linspace(1, length(cmap), noElements));
% 
% % Setup figure
% tiledlayout(1, 2);
% width = 6;
% height = 2;
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);
% 
% % Plotting times
% tPlots = [0.1, 0.2, 0.575];
% timesteps = cast(impact_timestep + tPlots / DELTA_T, "int32");
% 
% % Plot membrane in time
% for k = timesteps
%     t = ts(k);
% 
%     % Loop over cases
%     for paramIdx = noElements : -1 : 1
%         % Set param dir
%         param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", vary_dir, ...
%             ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));
% 
%         % Set line color
%         lineColor = cmap(DNSLineColorIdxs(paramIdx), :);
% 
%         % Set display name
%         displayName = append("$\gamma = $", GAMMA_strs(paramIdx));
% 
%         % Load dns files
%         dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1);
%         dns_mat = importdata(dns_filename);
%         xs_dns = dns_mat(:, 1);
%         [sort_xs, sort_idxs] = sort(xs_dns);
%         ws_dns = dns_mat(sort_idxs, 2);
%         w_ts_dns = dns_mat(sort_idxs, 3);
%         ps_dns = dns_mat(sort_idxs, 4);
%         
%         % Plot membrane position
%         nexttile(1);
%         h(paramIdx) = plot(sort_xs, -ws_dns, 'LineWidth', lineWidth, ...
%             'Color', lineColor, 'DisplayName', displayName);
%         hold on;
% 
%         % Plot membrane derivative
%         nexttile(2);
%         plot(sort_xs, -w_ts_dns, 'LineWidth', lineWidth, 'Color', lineColor);
%         hold on;
%     end
% 
%     % Set figure options
%     nexttile(1);
%     hold off;
%     xlim([0, 16]);
%     xticks(xTicks);
%     ylim([-0.2, 0.045]);
%     grid on;
%     box on;
%     xlabel("$x$");
%     ylabel("$-w(x, t)$");
% 
%     nexttile(2);
%     hold off;
%     xlim([0, 16]);
%     xticks(xTicks);
%     ylim([-0.45, 0.15])
%     grid on;
%     box on;
%     xlabel("$x$");
%     ylabel("$-w_t(x, t)$");
% 
%     % Set legend options
% %     lh = legend(h(1 : noElements), ...
% %         'NumColumns', 3);
% %     lh.Layout.Tile = 'North'; 
% 
%     % Pause
%     pause(0.01);
% 
%     % Output
%     figname = append("MembraneFigures/MembraneParamVary/MembraneGAMMAVary_", num2str(t));
%     exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
% end
% 
% % Plot pressure on substrate
% k = timesteps(1)
% 
% close(figure(1));
% figure(1);
% width = 6;
% height = 3;
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);
% 
% % Loop over cases
% for paramIdx = 1 : noElements
%     % Set param dir
%     param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", vary_dir, ...
%         ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));
% 
%     % Set line color
%     lineColor = cmap(DNSLineColorIdxs(paramIdx), :);
% 
%     % Set display name
%     displayName = append("$\gamma = $", GAMMA_strs(paramIdx));
% 
%     % Load dns files
%     dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1);
%     dns_mat = importdata(dns_filename);
%     xs_dns = dns_mat(:, 1);
%     [sort_xs, sort_idxs] = sort(xs_dns);
%     ps_dns = dns_mat(sort_idxs, 4);
%     
%     % Plot membrane position
%     h(paramIdx) = plot(sort_xs, ps_dns, 'LineWidth', lineWidth, ...
%         'Color', lineColor, 'Displayname', displayName);
%     hold on;
% 
% end
% 
% % Plot stationary pressure
%  dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", stat_dns_dir, k - 1);
% dns_mat = importdata(dns_filename);
% xs_dns = dns_mat(:, 1);
% [sort_xs, sort_idxs] = sort(xs_dns);
% ps_dns = dns_mat(sort_idxs, 4);
% h(noElements + 1) = plot(sort_xs, ps_dns, 'LineWidth', 1.25 * lineWidth, ...
%         'Color', 'black', 'Linestyle', ':','Displayname', 'Stationary');
% 
% 
% % Set figure options
% hold off;
% xlim([0, 1]);
% ylim('padded')
% grid on;
% box on;
% xlabel("$x$");
% ylabel("$p(x, t)$");
% 
% 
% % Set legend options
% lh = legend(h(1 : noElements + 1), ...
%     'NumColumns', 3, 'Location', 'Northoutside');
% 
% % Pause
% pause(0.01);
% 
% % Output
% figname = "MembraneFigures/MembraneParamVary/PressureGAMMAVary";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);


%% New Gamma varying
GAMMA_strs = ["1.0", "31.62", "1000.0", "31620.0", "1000000.0"];
ALPHA_strs = ["0.1375", "0.1375", "0.1375", "0.1375", "0.1375"];
BETA_strs = ["0.0", "0.0", "0.0", "0.0", "0.0"];
noElements = length(GAMMA_strs);
vary_dir = new_gamma_dir;

% Line colours
DNSLineColorIdxs = floor(linspace(1, length(cmap), noElements));

% Setup figure
tiledlayout(1, 2);
width = 6;
height = 2;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Plotting times
tPlots = [0.1, 0.2, 0.575];
timesteps = cast(impact_timestep + tPlots / DELTA_T, "int32");

% Plot membrane in time
for k = timesteps
    t = ts(k);

    % Loop over cases
    for paramIdx = noElements : -1 : 1
        % Set param dir
        param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", vary_dir, ...
            ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));

        % Set line color
        lineColor = cmap(DNSLineColorIdxs(paramIdx), :);

        % Set display name
        displayName = append("$\gamma = $", GAMMA_strs(paramIdx));

        % Load dns files
        dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1)
        dns_mat = importdata(dns_filename);
        xs_dns = dns_mat(:, 1);
        [sort_xs, sort_idxs] = sort(xs_dns);
        ws_dns = dns_mat(sort_idxs, 2);
        w_ts_dns = dns_mat(sort_idxs, 3);
        ps_dns = dns_mat(sort_idxs, 4);
        
        % Plot membrane position
        nexttile(1);
        h(paramIdx) = plot(sort_xs, -ws_dns, 'LineWidth', lineWidth, ...
            'Color', lineColor, 'DisplayName', displayName);
        hold on;

        % Plot membrane derivative
        nexttile(2);
        plot(sort_xs, -w_ts_dns, 'LineWidth', lineWidth, 'Color', lineColor);
        hold on;
    end

    % Set figure options
    nexttile(1);
    hold off;
    xlim([0, 16]);
    xticks(xTicks);
    ylim([-0.5, 0.1]);
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$-w(x, t)$");

    nexttile(2);
    hold off;
    xlim([0, 16]);
    xticks(xTicks);
    ylim([-0.9, 0.4])
    yticks([-0.8, -0.4, 0, 0.4])
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$-w_t(x, t)$");

    % Set legend options
%     lh = legend(h(1 : noElements), ...
%         'NumColumns', 3);
%     lh.Layout.Tile = 'North'; 

    % Pause
    pause(0.01);

    % Output
    figname = append("MembraneFigures/MembraneParamVary/MembraneGAMMAVary_", num2str(t));
    exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
end

% Plot pressure on substrate
k = timesteps(1)

close(figure(1));
figure(1);
width = 6;
height = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Loop over cases
for paramIdx = 1 : noElements
    % Set param dir
    param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", vary_dir, ...
        ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));

    % Set line color
    lineColor = cmap(DNSLineColorIdxs(paramIdx), :);

    % Set display name
    displayName = append("$\gamma = $", GAMMA_strs(paramIdx));

    % Load dns files
    dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1);
    dns_mat = importdata(dns_filename);
    xs_dns = dns_mat(:, 1);
    [sort_xs, sort_idxs] = sort(xs_dns);
    ps_dns = dns_mat(sort_idxs, 4);
    
    % Plot membrane position
    h(paramIdx) = plot(sort_xs, ps_dns, 'LineWidth', lineWidth, ...
        'Color', lineColor, 'Displayname', displayName);
    hold on;

end

% Plot stationary pressure
 dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", stat_dns_dir, k - 1);
dns_mat = importdata(dns_filename);
xs_dns = dns_mat(:, 1);
[sort_xs, sort_idxs] = sort(xs_dns);
ps_dns = dns_mat(sort_idxs, 4);
h(noElements + 1) = plot(sort_xs, ps_dns, 'LineWidth', 1.25 * lineWidth, ...
        'Color', 'black', 'Linestyle', ':','Displayname', 'Stationary');


% Set figure options
hold off;
xlim([0, 1]);
ylim('padded')
grid on;
box on;
xlabel("$x$");
ylabel("$p(x, t)$");


% Set legend options
lh = legend(h(1 : noElements + 1), ...
    'NumColumns', 3, 'Location', 'Northoutside');

% Pause
pause(0.01);

% Output
figname = "MembraneFigures/MembraneParamVary/PressureGAMMAVary";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Beta varying
% BETA_strs = ["0.1", "3.162", "100", "3162", "100000"];
% GAMMA_strs = ["668.0", "668.0", "668.0", "668.0", "668.0"];
% ALPHA_strs = ["1.1", "1.1", "1.1", "1.1", "1.1"];
% 
% noElements = length(GAMMA_strs);
% vary_dir = beta_dir;
% 
% % Line colours
% DNSLineColorIdxs = floor(linspace(1, length(cmap), noElements));
% 
% % Plotting times
% tPlots = [0.1, 0.2, 0.575];
% timesteps = cast(impact_timestep + tPlots / DELTA_T, "int32");
% 
% % Setup figure
% tiledlayout(1, 2);
% width = 6;
% height = 2;
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);
% 
% % Plot membrane in time
% for k = timesteps
%     t = ts(k);
% 
%     % Loop over cases
%     for paramIdx = noElements : -1 : 1
%         % Set param dir
%         if paramIdx == 1
%             param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", delta_dir, ...
%                 "1.1", "0", "668.0");
% 
%             % Set display name
%             displayName = "$\beta = 0$";
%         else
%             param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", vary_dir, ...
%                 ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));
% 
%             % Set display name
%             displayName = append("$\beta = $", BETA_strs(paramIdx));
%         end
% 
%         % Set line color
%         lineColor = cmap(DNSLineColorIdxs(paramIdx), :);
% 
%         % Load dns files
%         dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1);
%         dns_mat = importdata(dns_filename);
%         xs_dns = dns_mat(:, 1);
%         [sort_xs, sort_idxs] = sort(xs_dns);
%         ws_dns = dns_mat(sort_idxs, 2);
%         w_ts_dns = dns_mat(sort_idxs, 3);
%         ps_dns = dns_mat(sort_idxs, 4);
%         
%         % Plot membrane position
%         nexttile(1);
%         h(paramIdx) = plot(sort_xs, -ws_dns, 'LineWidth', ...
%             lineWidth, 'Color', lineColor, 'DisplayName', displayName);
%         hold on;
% 
%         % Plot membrane derivative
%         nexttile(2);
%         plot(sort_xs, -w_ts_dns, 'LineWidth', lineWidth, 'Color', lineColor);
%         hold on;
% 
%         % Plot membrane pressure
% %         nexttile(3);
% %         plot(sort_xs, ps_dns, 'LineWidth', lineWidth, 'Color', lineColor);
% %         hold on;
%     end
% 
%     % Set figure options
%     nexttile(1);
%     hold off;
%     xlim([0, 16]);
%     xticks(xTicks);
%     ylim([-0.07, 0.01]);
%     grid on;
%     box on;
%     xlabel("$x$");
%     ylabel("$-w(x, t)$");
% 
%     nexttile(2);
%     hold off;
%     xlim([0, 16]);
%     xticks(xTicks);
%     ylim([-0.15, 0.1]);
%     yticks([-0.15, -0.1, -0.05, 0, 0.05, 0.1]);
%     grid on;
%     box on;
%     xlabel("$x$");
%     ylabel("$-w_t(x, t)$");
% 
%     % Set legend options
% %     lh = legend(h(1 : noElements), ...
% %         'NumColumns', 3);
% %     lh.Layout.Tile = 'North'; 
% 
%     % Pause
%     pause(0.01);
% 
%     % Output
%     figname = append("MembraneFigures/MembraneParamVary/MembraneBETAVary_", num2str(t));
%     exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
% end
% 
% % Plot pressure on substrate
% k = timesteps(1)
% 
% close(figure(1));
% figure(1);
% width = 6;
% height = 3;
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);
% 
% % Loop over cases
% for paramIdx = 1 : noElements
%     % Set param dir
%     if paramIdx == 1
%         param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", delta_dir, ...
%             "1.1", "0", "668.0");
% 
%         % Set display name
%         displayName = "$\beta = 0$";
%     else
%         param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", vary_dir, ...
%             ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));
% 
%         % Set display name
%         displayName = append("$\beta = $", BETA_strs(paramIdx));
%     end
% 
%     % Set line color
%     lineColor = cmap(DNSLineColorIdxs(paramIdx), :);
% 
%     % Load dns files
%     dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1);
%     dns_mat = importdata(dns_filename);
%     xs_dns = dns_mat(:, 1);
%     [sort_xs, sort_idxs] = sort(xs_dns);
%     ps_dns = dns_mat(sort_idxs, 4);
%     
%     % Plot membrane position
%     h(paramIdx) = plot(sort_xs, ps_dns, 'LineWidth', lineWidth, ...
%         'Color', lineColor, 'Displayname', displayName);
%     hold on;
% 
% end
% 
% % Plot stationary pressure
%  dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", stat_dns_dir, k - 1);
% dns_mat = importdata(dns_filename);
% xs_dns = dns_mat(:, 1);
% [sort_xs, sort_idxs] = sort(xs_dns);
% ps_dns = dns_mat(sort_idxs, 4);
% h(noElements + 1) = plot(sort_xs, ps_dns, 'LineWidth', 1.25 * lineWidth, ...
%         'Color', 'black', 'Linestyle', ':','Displayname', 'Stationary');
% 
% 
% % Set figure options
% hold off;
% xlim([0, 1]);
% ylim('padded')
% grid on;
% box on;
% xlabel("$x$");
% ylabel("$p(x, t)$");
% 
% 
% % Set legend options
% lh = legend(h(1 : noElements + 1), ...
%     'NumColumns', 3, 'Location', 'Northoutside');
% 
% % Pause
% pause(0.01);
% 
% % Output
% figname = "MembraneFigures/MembraneParamVary/PressureBETAVary";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% NEW Beta varying
BETA_strs = ["0", "5.335", "94.87", "1687.0", "30000.0"];
GAMMA_strs = ["1.305", "1.305", "1.305", "1.305", "1.305"];
ALPHA_strs = ["0.1375", "0.1375", "0.1375", "0.1375", "0.1375"];

noElements = length(GAMMA_strs);
vary_dir = new_beta_dir;

% Line colours
DNSLineColorIdxs = floor(linspace(1, length(cmap), noElements));

% Plotting times
tPlots = [0.1, 0.2, 0.575];
timesteps = cast(impact_timestep + tPlots / DELTA_T, "int32");

% Setup figure
tiledlayout(1, 2);
width = 6;
height = 2;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Plot membrane in time
for k = timesteps
    t = ts(k);

    % Loop over cases
    for paramIdx = noElements : -1 : 1
        % Set param dir
        if paramIdx == 1
            param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", delta_dir, ...
                ALPHA_strs(paramIdx), "0", GAMMA_strs(paramIdx));

            % Set display name
            displayName = "$\beta = 0$";
        else
            param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", vary_dir, ...
                ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));

            % Set display name
            displayName = append("$\beta = $", BETA_strs(paramIdx));
        end

        % Set line color
        lineColor = cmap(DNSLineColorIdxs(paramIdx), :);

        % Load dns files
        dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1);
        dns_mat = importdata(dns_filename);
        xs_dns = dns_mat(:, 1);
        [sort_xs, sort_idxs] = sort(xs_dns);
        ws_dns = dns_mat(sort_idxs, 2);
        w_ts_dns = dns_mat(sort_idxs, 3);
        ps_dns = dns_mat(sort_idxs, 4);
        
        % Plot membrane position
        nexttile(1);
        h(paramIdx) = plot(sort_xs, -ws_dns, 'LineWidth', ...
            lineWidth, 'Color', lineColor, 'DisplayName', displayName);
        hold on;

        % Plot membrane derivative
        nexttile(2);
        plot(sort_xs, -w_ts_dns, 'LineWidth', lineWidth, 'Color', lineColor);
        hold on;

        % Plot membrane pressure
%         nexttile(3);
%         plot(sort_xs, ps_dns, 'LineWidth', lineWidth, 'Color', lineColor);
%         hold on;
    end

    % Set figure options
    nexttile(1);
    hold off;
    xlim([0, 16]);
    xticks(xTicks);
    ylim([-0.5, 0.1]);
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$-w(x, t)$");

    nexttile(2);
    hold off;
    xlim([0, 16]);
    xticks(xTicks);
    ylim([-0.9, 0.4]);
    yticks([-0.8, 0, 0.4]);
    grid on;
    box on;
    xlabel("$x$");
    ylabel("$-w_t(x, t)$");

    % Set legend options
%     lh = legend(h(1 : noElements), ...
%         'NumColumns', 3);
%     lh.Layout.Tile = 'North'; 

    % Pause
    pause(0.01);

    % Output
    figname = append("MembraneFigures/MembraneParamVary/MembraneBETAVary_", num2str(t));
    exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
end

% Plot pressure on substrate
k = timesteps(1)

close(figure(1));
figure(1);
width = 6;
height = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

% Loop over cases
for paramIdx = 1 : noElements
    % Set param dir
    if paramIdx == 1
        param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", delta_dir, ...
            ALPHA_strs(paramIdx), "0", GAMMA_strs(paramIdx));

        % Set display name
        displayName = "$\beta = 0$";
    else
        param_dir = sprintf("%s/alpha_%s-beta_%s-gamma_%s", vary_dir, ...
            ALPHA_strs(paramIdx), BETA_strs(paramIdx), GAMMA_strs(paramIdx));

        % Set display name
        displayName = append("$\beta = $", BETA_strs(paramIdx));
    end

    % Set line color
    lineColor = cmap(DNSLineColorIdxs(paramIdx), :);

    % Load dns files
    dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", param_dir, k - 1);
    dns_mat = importdata(dns_filename);
    xs_dns = dns_mat(:, 1);
    [sort_xs, sort_idxs] = sort(xs_dns);
    ps_dns = dns_mat(sort_idxs, 4);
    
    % Plot membrane position
    h(paramIdx) = plot(sort_xs, ps_dns, 'LineWidth', lineWidth, ...
        'Color', lineColor, 'Displayname', displayName);
    hold on;

end

% Plot stationary pressure
 dns_filename = sprintf("%s/membrane_outputs/membrane_arr_%d.txt", stat_dns_dir, k - 1);
dns_mat = importdata(dns_filename);
xs_dns = dns_mat(:, 1);
[sort_xs, sort_idxs] = sort(xs_dns);
ps_dns = dns_mat(sort_idxs, 4);
h(noElements + 1) = plot(sort_xs, ps_dns, 'LineWidth', 1.25 * lineWidth, ...
        'Color', 'black', 'Linestyle', ':','Displayname', 'Stationary');


% Set figure options
hold off;
xlim([0, 1]);
ylim('padded')
grid on;
box on;
xlabel("$x$");
ylabel("$p(x, t)$");


% Set legend options
lh = legend(h(1 : noElements + 1), ...
    'NumColumns', 3, 'Location', 'Northoutside');

% Pause
pause(0.01);

% Output
figname = "MembraneFigures/MembraneParamVary/PressureBETAVary";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
