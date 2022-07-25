%% force_validation.m
% Validates the force in the DNS by comparing it for multiple levels

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/Imposed/");
addpath("../Analytical_Scripts/Imposed/Forces");

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;



%% Figure options
fontsize = 12;
lineWidth = 1.25;
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize', fontsize);
set(0, 'DefaultTextFontSize', fontsize);
set(0,'defaultLegendFontSize', fontsize, 'DefaultLegendFontSizeMode','manual');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(0, 'DefaultLegendInterpreter', 'latex');

%% Load in analytical parameters
% analyticalType = "flatDNS";
% [epsilon, q, omega, p, L] = substrateparameters(analyticalType);

%% Parameters
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
tsAnalytical = 1e-9 : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);
freq = 25; % Frequency

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
types = ["stationary_plate", "imposed_plate_0.00125"];
levels = 10 : 14;

%% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels - 1, 3);

for q = 1 : no_levels - 1
    colors(q, :) = cmap(color_idxs(q), :);
end

%% Norms matrix


%% Plot all cases
for axi = [0, 1]
% for axi = 0
    tileFig = tiledlayout(1, 2, 'Tilespacing', 'compact', 'padding', 'compact');
    
    % Norms matrix
    norms = zeros(length(levels) - 1, 2);
    for typeIdx = 1 : length(types)
        type = types(typeIdx);
        
        tileFig;
        % Save for parent directory
        parent_dir = sprintf("%s/%s_maxlevel_validation/axi_%d", ...
            master_dir,type, axi);
        
        %% Analytical solutions
        % Set analytical type
        if type == "stationary_plate"
            analyticalType = "stationary";
            titleStr = "(a) Stationary substrate setup.";
        else
            analyticalType = "flatDNS";
            titleStr = "(b) Moving frame setup.";
        end
        
        % Set analytical dimension
        if axi == 0
            dimension = "2D";
        else
            dimension = "axi";
        end
        
        % Load analytical substrate functions and parameters
        [epsilon, ~, ~, ~, ~] = substrateparameters(analyticalType);
        SubstrateFunctions = substratefunctions(analyticalType, dimension);
        
        % Find analytical times
        tsAnalytical = (1e-9 : DELTA_T : T_MAX - IMPACT_TIME) / epsilon^2;
        
        % Load analytical turnover points
        dsAnalytical = SubstrateFunctions.d(tsAnalytical);
        
        % Find times restricted to where where turnover point less that 1
        tsAnalyticalRestricted = tsAnalytical(epsilon * dsAnalytical <= 1);
        
        % Find analytical force
        [Fs_composite, Fs_outer, ~] ...
            = substrateforce(tsAnalyticalRestricted, SubstrateFunctions);
        
        %% Plot DNS solutions
        nexttile(typeIdx);
        hold on;
        
        % Matrix to hold the turnover points
        dns_forces = zeros(NO_TIMESTEPS, length(levels));
        
        % Loop over the levels
        for level_idx = 1 : length(levels)
            level = levels(level_idx);
            level_dir = sprintf("%s/max_level_%d", parent_dir, level);

            forces_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", level_dir));
            ts = forces_mat(1 : NO_TIMESTEPS, 1) - IMPACT_TIME;
            Fs = forces_mat(1 : NO_TIMESTEPS, 3);

            dns_forces(:, level_idx) = Fs;

%             scatter(ts(1 : freq : end), Fs(1 : freq : end), [], cmap(color_idxs(level_idx), :), ...
%                 'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
            h(level_idx) = plot(ts, Fs, 'color', cmap(color_idxs(level_idx), :), ...
                'linewidth', lineWidth, 'Displayname', sprintf("$m$ = %d", level));
        end
        
        %% Determine norms
        Fs_max = dns_forces(:, length(levels));

        for level_idx = 1 : length(levels) - 1
            level = levels(level_idx);

            % Determines L2-norm difference
            diff = dns_forces(:, level_idx) - Fs_max;
            norms(level_idx, typeIdx) = sqrt(sum(diff.^2) / length(diff));
        end
        
        %% Plot analytical solutions
        % Plot empty white line
        h(length(levels) + 1) = scatter(10, 10, [], 'white', 'Displayname', '');
        
        % Outer force
        h(length(levels) + 2) = plot(tsAnalyticalRestricted * epsilon^2, Fs_outer, ...
            'linestyle', '--', 'linewidth', lineWidth, ...
            'color', 'black', 'displayname', 'Analytical (leading-order)');
        
        % Composite force
        h(length(levels) + 3) = plot(tsAnalyticalRestricted * epsilon^2, Fs_composite, ...
            'linewidth', lineWidth, ...
            'color', 'black', 'displayname', 'Analytical (composite)');
        
        % Plot vertical line where analytical solution ends
        xline(tsAnalyticalRestricted(end) * epsilon^2, 'linestyle', '--', 'color', 0.5 * [1 1 1]);
        
        grid on;
        if axi == 1
            ylim([-0.5, 7]);
        else
            ylim([-0.5, 8]);
        end
        xlim([-0.2, 0.8]);
        xticks(-0.4 : 0.2 : 0.8);
        
        xlabel("$t$");
        ylabel("$F_m(t)$");
        
        %% Create title
        title(titleStr, 'Fontsize', fontsize);
        set(gca, 'TitleFontSizeMultiplier', 1);
    end
   
    
    %% Set figure size
    lh = legend(h(1 : length(levels) + 3), 'interpreter', 'latex', 'Numcolumns', 3);
    lh.Layout.Tile = 'South'; 
    
    % Pixel options
%     x0=400;
%     y0=400;
%     height=700;
%     width=1150;
%     set(gcf,'position',[x0,y0,width,height]);
    
    % Inches options
    width = 6;
    height = 3.65;
    set(gcf,'units', 'inches', 'position',[0.5 * width, 0.5 * height, width, height]);
    
    %% Plot L2 norm error
    for typeIdx = 1 : length(types)
        
        insetWidth = 0.10;
        insetHeight = 0.16;
        if axi == 1
            verticalPos = 0.7;
            if typeIdx == 1
                horizontalPos = 0.362;
            else
                horizontalPos = 0.848;
            end
        else
            verticalPos = 0.7;
            if typeIdx == 1
                horizontalPos = 0.362;
            else
                horizontalPos = 0.848;
            end
        end
        
        axes('Position',[horizontalPos, verticalPos, insetWidth, insetHeight]);

        box on
        hold on;

        plot(levels(1 : end - 1), norms(:, typeIdx), 'color', 'black', ...
            'linewidth', lineWidth);

        sz = 35;
        scatter(levels(1 : end - 1), norms(:, typeIdx), sz, colors, 'filled');

        set(gca, 'yscale', 'log');
        xlim([min(levels) - 0.5, max(levels) - 0.5])
        ylim([10^-2, 1]);
        xticks(levels(1 : end - 1));

        grid on;
        xlabel("$m$");
        ylabel("$||F_m - F_{14}||_2$");
        hold off;
    end
    
    %% Save figure
    set(gcf, 'Renderer', 'Painters');
    pause(1.5);

    if axi == 0
        axiStr = "2D";
    else
        axiStr = "Axi";
    end
    figname = sprintf("dns_validation_figures/forces/DNSForce_%s", axiStr);
        
    exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
    savefig(gcf, sprintf("%s.fig", figname));
end