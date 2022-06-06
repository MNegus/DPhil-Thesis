%% force_validation.m
% Validates the force in the DNS by comparing it for multiple levels

clear;
close all;

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

fontsize = 22;

%% Parameters
L = 2;
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);
freq = 25; % Frequency
a = 0.05; % Acceleration parameter

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

%% Plot all cases

for type = types
    
    % Create parent_dirs
    coeffs = [2, sqrt(3)]; % Coeff of ds_analytical

    parent_dirs = strings(2, 1);
    for axi = [0, 1]
        parent_dirs(axi + 1) = sprintf("%s/%s_maxlevel_validation/axi_%d", ...
            master_dir,type, axi);
    end

    for parentIdx = 1 : length(parent_dirs)
        % Save for parent directory
        parent_dir = parent_dirs(parentIdx);
        axi = parentIdx - 1;
        
        % Analytical soltution
        if type == "stationary"
            ds_analytical = coeffs(parentIdx) * sqrt(ts_analytical);
        else
            ds_analytical = coeffs(parentIdx) * sqrt(ts_analytical ...
            - 0.5 * a * ts_analytical.^2);
        end
        
        % Create figure
        close(figure(1));
        figure(1);
        hold on;
        
        % Matrix to hold the turnover points
        dns_forces = zeros(NO_TIMESTEPS, length(levels));
        
        % Loop over the levels
        for level_idx = 1 : length(levels)
            level = levels(level_idx);
            level_dir = sprintf("%s/max_level_%d", parent_dir, level);

            level
            forces_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", level_dir));
            ts = forces_mat(1 : NO_TIMESTEPS, 1) - IMPACT_TIME;
            Fs = forces_mat(1 : NO_TIMESTEPS, 3);

            dns_forces(:, level_idx) = Fs;

%             scatter(ts(1 : freq : end), Fs(1 : freq : end), [], cmap(color_idxs(level_idx), :), ...
%                 'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
            plot(ts, Fs, 'color', cmap(color_idxs(level_idx), :), ...
                'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
        end
        
        % Analytical solution
%         plot(ts_analytical, ds_analytical, 'linestyle', '--', 'linewidth', 2, ...
%             'color', 'black', 'displayname', 'Analytical');
        
        grid on;
        if axi == 1
            ylim([-0.4, 3.5]);
        else
            ylim([-0.8, 7]);
        end
        xlim([-0.3, 0.8]);
%         xticks(-0.125 : 0.125 : 0.8);
        
        tix=get(gca,'ytick')';
        set(gca,'yticklabel',num2str(tix,'%.1f'));
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
        xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
        ylabel("$F_m(t)$", 'interpreter', 'latex');
        
        
        %% Plot L2 norm error
        if axi == 1
            legend('location', 'northwest', 'interpreter', 'latex');
            axes('Position',[.58 .24 .3 .3]);
        else
            legend('location', 'southeast', 'interpreter', 'latex');
            axes('Position',[.58 .60 .3 .3]);
        end
        box on
        hold on;
        Fs_max = dns_forces(:, length(levels));


        norms = zeros(length(levels) - 1, 1);
        for level_idx = 1 : length(levels) - 1
            level = levels(level_idx);

            % Determines L2-norm difference
            diff = dns_forces(:, level_idx) - Fs_max;
            norms(level_idx) = sqrt(sum(diff.^2) / length(diff));
        end
        plot(levels(1 : end - 1), norms, 'color', 'black', 'linewidth', 2);

        sz = 100;
        scatter(levels(1 : end - 1), norms, sz, colors, 'filled');

        set(gca, 'yscale', 'log');
        xlim([min(levels) - 0.5, max(levels) - 0.5])
        ylim([10^-2, 1]);
        xticks(levels(1 : end - 1));
%         yticks([10^-3, 10^-2, 10^-1]);

        grid on;
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
        xlabel("$m$", 'interpreter', 'latex');
        ylabel("$||F_m - F_{14}||_2$", 'interpreter', 'latex');

        %% Overal figure settings
        x0=400;
        y0=400;
        height=800;
        width=600;

        set(gcf,'position',[x0,y0,width,height]);
        set(gcf, 'Renderer', 'Painters');
        pause(1.5);

        figname = sprintf("dns_validation_figures/forces/DNSForce_%s_axi_%d", type, axi);
        
        exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
        savefig(gcf, sprintf("%s.fig", figname));
        
    end
end