%% bubblet_validation.m
% Validates the entrapped bubble size in the DNS by comparing it for 
% multiple levels

clear;
close all;

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

fontsize = 22;

%% Parameters
L = 3;
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);
freq = 100; % Frequency

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
        % Save for parent directory
        parent_dir = parent_dirs(parentIdx);
        
        % Create figure
        close(figure(1));
        figure(1);
        hold on;
        
        % Matrix to hold the turnover points
        dns_areas = zeros(NO_TIMESTEPS, length(levels));
        
        % Loop over the levels
        for level_idx = 1 : length(levels)
            level = levels(level_idx);
            level_dir = sprintf("%s/max_level_%d", parent_dir, level);

            
            if type == "membrane"
                bubbles_mat = readmatrix(sprintf("%s/output.txt", level_dir));
                ts = bubbles_mat(1 : NO_TIMESTEPS, 1) - IMPACT_TIME;
                areas = bubbles_mat(1 : NO_TIMESTEPS, 3);
            else
                bubbles_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", level_dir));
                ts = bubbles_mat(1 : NO_TIMESTEPS, 1) - IMPACT_TIME;
                areas = bubbles_mat(1 : NO_TIMESTEPS, 9);
                
                if parentIdx == 1
                    areas = areas / (2 * pi); % 2D case
                else
                    areas = areas / (4 * pi / 3);
                end
            end

            largeIdxs = areas > 1;
            areas(largeIdxs) = zeros(sum(largeIdxs), 1);
            
            dns_areas(:, level_idx) = areas;

%             scatter(ts(1 : freq : end), Fs(1 : freq : end), [], cmap(color_idxs(level_idx), :), ...
%                 'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
            scatter(ts(1 : freq : end), areas(1 : freq : end), [],  cmap(color_idxs(level_idx), :), ...
                'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
        end
        
        grid on;
%         ylim([0, 1.6]);
        xlim([-0.2, 0.7]);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
        xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
        ylabel("$A_m(t)$", 'interpreter', 'latex');
        legend('location', 'northwest', 'interpreter', 'latex');
        
        %% Plot L2 norm error
        axes('Position',[.34 .6 .3 .3])
        box on
        hold on;
        ds_max = dns_areas(:, length(levels));


        norms = zeros(length(levels) - 1, 1);
        for level_idx = 1 : length(levels) - 1
            level = levels(level_idx);

            % Determines L2-norm difference
            diff = dns_areas(:, level_idx) - ds_max;
            norms(level_idx) = sqrt(sum(diff.^2) / length(diff));
        end
        plot(levels(1 : end - 1), norms, 'color', 'black', 'linewidth', 2);

        sz = 100;
        scatter(levels(1 : end - 1), norms, sz, colors, 'filled');

        set(gca, 'yscale', 'log');
        xlim([min(levels) - 0.5, max(levels) - 0.5])
%         ylim([10^-3, 10^(-0.5)]);
%         xticks(levels(1 : end - 1));
%         yticks([10^-3, 10^-2, 10^-1]);

        grid on;
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
        xlabel("$m$", 'interpreter', 'latex');
        ylabel("$||A_m - d_{14}||_2$", 'interpreter', 'latex');

        %% Overal figure settings
        x0=400;
        y0=400;
        height=800;
        width=600;

        set(gcf,'position',[x0,y0,width,height]);
        set(gcf, 'Renderer', 'Painters');
        pause(1.5);

        if parentIdx == 1
            figname = "dns_validation_figures/DNSBubbles_2D";
        else
            figname = "dns_validation_figures/DNSBubbles_axi";
        end
%         figname = sprintf("dns_validation_figures/bubbles/DNSBubbles_%s_idx_%d", type, parentIdx);
        exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
        savefig(gcf, sprintf("%s.fig", figname));
        
    end
end
    