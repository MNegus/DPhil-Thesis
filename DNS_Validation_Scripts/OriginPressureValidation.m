%% OriginPressureValidation.m
% Validates the force in the DNS by comparing it for multiple levels

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/Imposed/");
addpath("../Analytical_Scripts/Imposed/Pressures");

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
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);
freq = 100; % Frequency of scatter

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
type = ["imposed_0.05_omega_4"];
levels = 10 : 14;
imposed_coeffs = ["0", "0.05"]; % Coefficients of imposed substrate

%% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels - 1, 3);

for q = 1 : no_levels - 1
    colors(q, :) = cmap(color_idxs(q), :);
end

%% Plot all cases
% Set up tiled figure
tileFig = tiledlayout(1, 2, 'Tilespacing', 'compact', 'padding', 'compact');

% Norms matrix
norms = zeros(length(levels) - 1, 2);

% Loop over types (stationary or moving) 
for imposedIdx = 1 : 2
    imposed_coeff = imposed_coeffs(imposedIdx); % Save coefficient

    tileFig; % Change to tiled figure

    % Save parent directory
    parent_dir = sprintf("%s/%s_maxlevel_validation/imposed_coeff_%s", ...
        master_dir, type, imposed_coeff);

    %% Determine analytical solutions
    % Set analytical type and title string
    if imposed_coeff == "0"
        analyticalType = "stationary";
        titleStr = "(a) Stationary substrate setup.";
    else
        analyticalType = "curvedDNS";
        titleStr = "(b) Linearised boundary setup.";
    end
    
    % Load analytical substrate functions and parameters
    [epsilon, ~, ~, ~, ~] = substrateparameters(analyticalType);
    SubstrateFunctions = substratefunctions(analyticalType, "2D");

    % Find analytical times
    tsAnalytical = (1e-9 : DELTA_T : T_MAX - IMPACT_TIME) / epsilon^2;

    % Load analytical turnover points
    dsAnalytical = SubstrateFunctions.d(tsAnalytical);

    % Find times restricted to where where turnover point less that 1
    tsAnalyticalRestricted = tsAnalytical(epsilon * dsAnalytical <= 1);

    % Determine pressure at the origin
    p0sAnalytical = zeros(size(tsAnalyticalRestricted));
    for k = 1 : length(tsAnalyticalRestricted)
        t = tsAnalyticalRestricted(k);
        p0sAnalytical(k) = outerpressure(0, t, SubstrateFunctions);
    end
    
    % Save analytical solutions in structure array
    if imposed_coeff == "0"
        analyticalStruct.tsStationary = tsAnalyticalRestricted;
        analyticalStruct.p0sStationary = p0sAnalytical;
    else
        analyticalStruct.tsMoving = tsAnalyticalRestricted;
        analyticalStruct.p0sMoving = p0sAnalytical;
    end

    %% Plot DNS solutions
    % Move to relevant tile
    nexttile(imposedIdx);
    hold on;

    % Matrix to hold the force locations
    dns_pressures = zeros(NO_TIMESTEPS, length(levels));

    % Loop over the levels
    for level_idx = 1 : length(levels)
        % Save level and directory
        level = levels(level_idx);
        level_dir = sprintf("%s/max_level_%d", parent_dir, level);

        % Load matrix containing forces
        pressures_mat = readmatrix(sprintf("%s/output.txt", level_dir));
        ts = pressures_mat(1 : NO_TIMESTEPS, 1) - IMPACT_TIME;
        p0s = pressures_mat(1 : NO_TIMESTEPS, 3);

        dns_pressures(:, level_idx) = p0s; % Save forces into matrix

        % Plot DNS pressures
        sz = 15;
        h(level_idx) = scatter(ts(1 : freq : end), p0s(1 : freq : end), ...
            sz, cmap(color_idxs(level_idx), :), 'linewidth', lineWidth, ...
            'Displayname', sprintf("$m$ = %d", level));
    end

    %% Determine norms
    p0s_max = dns_pressures(:, length(levels)); % Force of maximum level

    % Loop over levels
    for level_idx = 1 : length(levels) - 1
        level = levels(level_idx);

        % Determines L2-norm difference
        diff = dns_pressures(:, level_idx) - p0s_max;
        norms(level_idx, imposedIdx) = sqrt(sum(diff.^2) / length(diff));
    end

    %% Plot analytical solution
    h(length(levels) + 1) = plot(tsAnalyticalRestricted * epsilon^2, p0sAnalytical, ...
        'linestyle', '--', 'linewidth', lineWidth, ...
        'color', 'black', 'displayname', 'Analytical');

    % Plot vertical line where analytical solution ends
    xline(tsAnalyticalRestricted(end) * epsilon^2, ...
        'linestyle', '--', 'color', 0.5 * [1 1 1]);

    %% Tile figure options
    grid on;
    set(gca, 'YScale', 'log')
    
    % Set axes limits and ticks
    xlim([-0.2, 0.8]);
    xticks(-0.4 : 0.2 : 0.8);
    ylim([0, 100]);

    % Axes labels
    xlabel("$t$");
    ylabel("$p_m(0, t)$");

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
height = 4.5;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

%% Plot L2 norm error as inset plot
% Loop over types
for imposedIdx = 1 : 2

    % Define axes options for the inset
    insetWidth = 0.11;
    insetHeight = 0.16;
    verticalPos = 0.75;
    if imposedIdx == 1
        horizontalPos = 0.35;
    else
        horizontalPos = 0.835;
    end

    % Define inset axes
    axes('Position',[horizontalPos, verticalPos, insetWidth, insetHeight]);
    box on;
    grid on;
    hold on;

    % Plot black line for norms
    plot(levels(1 : end - 1), norms(:, imposedIdx), 'color', 'black', ...
        'linewidth', lineWidth);

    % Plot scatter for norms solution
    sz = 35;
    scatter(levels(1 : end - 1), norms(:, imposedIdx), sz, colors, 'filled');
    hold off;

    % Set axes options
    set(gca, 'yscale', 'log');
    xlim([min(levels) - 0.5, max(levels) - 0.5])
    ylim([0.1, 2]);
    
    xticks(levels(1 : end - 1));
%     yticks([10^-3, 10^-2, 10^-1]);

    % Set axes labels
    xlabel("$m$");
    ylabel("$||p_m - p_{14}||_2$");
end

%% Save figure
set(gcf, 'Renderer', 'Painters');
pause(1.5);

% Figure name
figname = "dns_validation_figures/pressures/DNSOriginPressure";

% Export figure
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
savefig(gcf, sprintf("%s.fig", figname));
