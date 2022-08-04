%% ForceValidation.m
% Validates the turnover point in the DNS by comparing it for multiple 
% levels.

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/Imposed/");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;

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
DELTA_T = 10 * 1e-4;
IMPACT_TIME = 0.125;
freq = 7; % Frequency to plot scatters

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
type = "stationary_plate";
levels = 10 : 14;

%% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels - 1, 3);

for q = 1 : no_levels - 1
    colors(q, :) = cmap(color_idxs(q), :);
end

%% Plot all cases
tileFig = tiledlayout(1, 2, 'Tilespacing', 'compact', 'padding', 'compact');

% Norms matrix
norms = zeros(length(levels) - 1, 2);

% Loop over dimensions
for axi = [0, 1]

    tileFig; % Change to tiled figure

    % Save parent directory
    parent_dir = sprintf("%s/%s_maxlevel_validation/axi_%d", ...
        master_dir, type, axi);

    %% Determine analytical solutions

    % Set dimension specific parameters
    analyticalType = "stationary";
    if axi == 0
        T_MAX = 0.45;
        dimension = "2D";
        titleStr = "(a) Two-dimensional case.";
    else
        T_MAX = 0.55;
        dimension = "axi";
        titleStr = "(b) Axisymmetric case.";
    end
    
    % Time variables
    ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
    MAX_TIMESTEP = T_MAX / DELTA_T;
    NO_TIMESTEPS = length(ts);

    % Load analytical substrate functions and parameters
    [epsilon, ~, ~, ~, ~] = substrateparameters(analyticalType);
    SubstrateFunctions = substratefunctions(analyticalType, dimension);

    % Find analytical times
    tsAnalytical = (1e-9 : DELTA_T : T_MAX - IMPACT_TIME) / epsilon^2;

    % Load analytical turnover points
    dsAnalytical = SubstrateFunctions.d(tsAnalytical);

    % Find times restricted to where where turnover point less that 1
    tsAnalyticalRestricted = tsAnalytical(epsilon * dsAnalytical <= 1);
    dsAnalyticalRestricted = dsAnalytical(epsilon * dsAnalytical <= 1);

    %% Plot DNS solutions
    % Move to relevant tile
    nexttile(axi + 1);
    hold on;

    % Matrix to hold the force locations
    dns_turnovers = zeros(NO_TIMESTEPS, length(levels));

    % Loop over the levels
    for level_idx = 1 : length(levels)
        % Save level and directory
        level = levels(level_idx);
        level_dir = sprintf("%s/max_level_%d", parent_dir, level);

        % Load matrix containing turnover points
        turnover_mat = readmatrix(sprintf("%s/turnover_points.csv", level_dir));
        ds = turnover_mat(1 : NO_TIMESTEPS, 2);

        dns_turnovers(:, level_idx) = ds; % Save turnover points into matrix

        % Plot DNS turnover points
        sz = 15;
        h(level_idx) = scatter(ts(1 : freq : end), ds(1 : freq : end), ...
            sz, cmap(color_idxs(level_idx), :), 'linewidth', lineWidth, ...
            'Displayname', sprintf("$m$ = %d", level));
    end

    %% Determine norms
    ds_max = dns_turnovers(:, length(levels)); % Force of maximum level

    % Loop over levels
    for level_idx = 1 : length(levels) - 1
        level = levels(level_idx);

        % Determines L2-norm difference
        diff = dns_turnovers(:, level_idx) - ds_max;
        norms(level_idx, axi + 1) = sqrt(sum(diff.^2) / length(diff));
    end

    %% Plot analytical solutions
    % Plot outer force
    h(length(levels) + 1) = plot(tsAnalyticalRestricted * epsilon^2, ...
        epsilon * dsAnalyticalRestricted, ...
        'linestyle', '--', 'linewidth', lineWidth, ...
        'color', 'black', 'displayname', 'Analytical');

    % Plot vertical line where analytical solution ends
    xline(tsAnalyticalRestricted(end) * epsilon^2, ...
        'linestyle', '--', 'color', 0.5 * [1 1 1]);

    %% Tile figure options
    grid on;

    % Set axes limits and ticks
    xlim([-0.2, 0.45]);
    xticks(-0.4 : 0.2 : 0.8);
    ylim([0, 1.2]);

    % Axes labels
    xlabel("$t$");
    ylabel("$d_m(t)$");

    % Create title
    title(titleStr, 'Fontsize', fontsize);
    set(gca, 'TitleFontSizeMultiplier', 1);
end

%% Set options for entire figure
% Set legend
lh = legend(h(1 : length(levels) + 1), 'interpreter', 'latex', 'Numcolumns', 3);
lh.Layout.Tile = 'South'; 

% Set size in inches
width = 6;
height = 3.65;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

%% Plot L2 norm error as inset plot
% Loop over types
for axi = [0, 1]

    % Define axes options for the inset
    insetWidth = 0.11;
    insetHeight = 0.16;
    verticalPos = 0.72;
    if axi == 0
        horizontalPos = 0.17;
    else
        horizontalPos = 0.65;
    end

    % Define inset axes
    axes('Position',[horizontalPos, verticalPos, insetWidth, insetHeight]);
    box on;
    grid on;
    hold on;

    % Plot black line for norms
    plot(levels(1 : end - 1), norms(:, axi + 1), 'color', 'black', ...
        'linewidth', lineWidth);

    % Plot scatter for norms solution
    sz = 35;
    scatter(levels(1 : end - 1), norms(:, axi + 1), sz, colors, 'filled');
    hold off;

    % Set axes options
    set(gca, 'yscale', 'log');
    xlim([min(levels) - 0.5, max(levels) - 0.5])
    ylim([10^-3, 10^(-0.5)]);
    
    xticks(levels(1 : end - 1));
    yticks([10^-3, 10^-2, 10^-1]);

    % Set axes labels
    xlabel("$m$");
    ylabel("$||d_m - d_{14}||_2$");
end

%% Save figure
set(gcf, 'Renderer', 'Painters');
pause(1.5);

% Figure name
figname = "dns_validation_figures/turnovers/DNSTurnover";

% Export figure
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
savefig(gcf, sprintf("%s.fig", figname));
