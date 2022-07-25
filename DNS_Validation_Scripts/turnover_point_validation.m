%% turnover_point_validation.m
% Validates the turnover point location in the DNS by comparing its value 
% for multiple levels

clear;
close all;

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

fontsize = 22;

%% Parameters
L = 3;
DELTA_T = 1e-3;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
% ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
% NO_TIMESTEPS = length(ts);
freq = 7; % Frequency
a = 0.05; % Acceleration parameter

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
types = ["stationary_plate"];
levels = 10 : 14;

%% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels - 1, 3);

for q = 1 : no_levels - 1
    colors(q, :) = cmap(color_idxs(q), :);
end

%% Plot all cases

for type = types
    
    % Create parent_dirs
    if (type == "stationary_plate" || type == "imposed_0.00125")
        coeffs = [2, sqrt(3)]; % Coeff of ds_analytical
        
        parent_dirs = strings(2, 1);
        for axi = [0, 1]
            parent_dirs(axi + 1) = sprintf("%s/%s_maxlevel_validation/axi_%d", ...
                master_dir,type, axi);
        end
    else
        parent_dirs = sprintf("%s/%s_maxlevel_validation/membrane_acc_0.05", ...
            master_dir, type); 
    end
    
    for parentIdx = 1 : length(parent_dirs)
%     for parentIdx = 1
        % Save for parent directory
        parent_dir = parent_dirs(parentIdx)
        
        % Analytical soltution
        if type == "membrane"
            ds_guess = 2 * sqrt(ts_analytical);
            k = pi / (2 * L);
            zero_fun = @(d) ts_analytical - d.^2 / 4 ...
                - 0.5 * a * ts_analytical.^2 .* besselj(0, k * d);
            ds_analytical = fsolve(zero_fun, ds_guess);
        elseif type == "stationary_plate"
            ds_analytical =  @(t) coeffs(parentIdx) * sqrt(t);
        else
            ds_analytical = @(t) coeffs(parentIdx) * sqrt(t);
        end
        
        % Create figure
        close(figure(1));
        figure(1);
        hold on;
        
        % Find the timestep where the maxlevel solution exceed d = 1
        level_dir = sprintf("%s/max_level_%d", parent_dir, max(levels));
        turnover_mat = readmatrix(sprintf("%s/turnover_points.csv", level_dir));
        ds = turnover_mat(:, 2);
        max_timestep = sum(ds <= 1) + 6 * freq;
        
        % Matrix to hold the turnover points
        dns_turnovers = zeros(max_timestep, length(levels));
        
        % Loop over the levels
        for level_idx = 1 : length(levels)
            level = levels(level_idx);
            level_dir = sprintf("%s/max_level_%d", parent_dir, level);

            turnover_mat = readmatrix(sprintf("%s/turnover_points.csv", level_dir));
            ts = turnover_mat(1 : max_timestep, 1) - IMPACT_TIME;
            size(ts)
            ds = turnover_mat(1 : max_timestep, 2);
            dns_turnovers(:, level_idx) = ds;

            scatter(ts(1 : freq : end), ds(1 : freq : end), [], cmap(color_idxs(level_idx), :), ...
                'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
        end
        
        % Analytical solution
        T_MAX = DELTA_T * max_timestep;
        ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
        plot(ts_analytical, ds_analytical(ts_analytical), 'linestyle', '--', 'linewidth', 2, ...
            'color', 'black', 'displayname', 'Analytical');
        
        grid on;
        xlim([-0.2, 0.45]);
        ylim([0, 1.7]);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
        xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
        ylabel("$d_m(t)$", 'interpreter', 'latex');
        legend('location', 'southeast', 'interpreter', 'latex');
        
        impact_idx = IMPACT_TIME / DELTA_T;
        plot(ts(impact_idx : end), ds(impact_idx : end));
        
        %% Plot L2 norm error
        axes('Position',[.34 .6 .3 .3])
        box on
        hold on;
        
        % Turnovers for max level
        ds_max = dns_turnovers(:, length(levels));
        

        norms = zeros(length(levels) - 1, 1);
        for level_idx = 1 : length(levels) - 1
            level = levels(level_idx);

            % Determines L2-norm difference
            diff = dns_turnovers(:, level_idx) - ds_max;
            norms(level_idx) = sqrt(sum(diff.^2) / length(diff));
        end
        plot(levels(1 : end - 1), norms, 'color', 'black', 'linewidth', 2);

        sz = 100;
        scatter(levels(1 : end - 1), norms, sz, colors, 'filled');

        set(gca, 'yscale', 'log');
        xlim([min(levels) - 0.5, max(levels) - 0.5])
        ylim([10^-3, 10^(-0.5)]);
        xticks(levels(1 : end - 1));
        yticks([10^-3, 10^-2, 10^-1]);

        grid on;
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
        xlabel("$m$", 'interpreter', 'latex');
        ylabel("$||d_m - d_{14}||_2$", 'interpreter', 'latex');

        %% Overal figure settings
        x0=400;
        y0=400;
        height=800;
        width=600;

        set(gcf,'position',[x0,y0,width,height]);
        pause(1.5);

%         figname = sprintf("dns_validation_figures/turnovers/DNSTurnover_%s_idx_%d", type, parentIdx);
%         set(gcf, 'Renderer', 'Painters');
%         exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
%         savefig(gcf, sprintf("%s.fig", figname));
        
    end
end
    