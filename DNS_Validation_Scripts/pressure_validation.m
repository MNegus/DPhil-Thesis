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
L = 2;
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.475;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);
freq = 5; % Frequency

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
% type = ["imposed_membrane_0.0025"];
type = ["EasyCase"];
levels = 10 : 13;

%% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels - 1, 3);

for q = 1 : no_levels - 1
    colors(q, :) = cmap(color_idxs(q), :);
end

%% Save analytical solution
epsilon = 1; % Analytical parameter for small time 
% a = 0.0025; % Imposed parameter
% [ws, w_ts, w_tts] = ImposedPlate(ts_analytical, a);
% lambda = pi / (2 * L);
% p0sAnalytical = AnalyticalPressure2D(ts_analytical, ws, w_ts, w_tts, lambda, epsilon);


%% Plot all cases
imposed_coeffs = ["0"; "0.00125"; "0.0025"];
for imposedIdx = 1:3
    
    imposed_coeff = imposed_coeffs(imposedIdx);
    
    
    
    if imposed_coeff == "0"
        parent_dir = "/scratch/negus/DNS_Chapter_validation/imposed_coeff_0";
        p0sAnalytical = 1 ./ sqrt(ts_analytical); 
    else
        parent_dir = sprintf("%s/%s_maxlevel_validation/imposed_coeff_%s", ...
        master_dir, type, imposed_coeff);
%         p0sAnalytical = AnalyticalPressure2D(ts_analytical, ws, w_ts, w_tts, lambda, epsilon);
        p0sAnalytical = 1 ./ sqrt(ts_analytical); 
    end

    % Create figure
    close(figure(imposedIdx));
    figure(imposedIdx);
    hold on;

    % Matrix to hold the turnover points
    dns_ps = zeros(NO_TIMESTEPS, length(levels));

    % Loop over the levels
    for level_idx = 1 : length(levels)
        level = levels(level_idx);
        level_dir = sprintf("%s/max_level_%d", parent_dir, level);

        pressures_mat = readmatrix(sprintf("%s/output.txt", level_dir));
        ts = pressures_mat(1 : NO_TIMESTEPS, 1) - IMPACT_TIME;
        p0s = pressures_mat(1 : NO_TIMESTEPS, 3);

        dns_ps(:, level_idx) = p0s;

        scatter(ts(1 : freq : end), p0s(1 : freq : end), [], cmap(color_idxs(level_idx), :), ...
            'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
        drawnow;
    end
    % Analytical solution
    plot(ts_analytical, p0sAnalytical, 'linestyle', '--', 'linewidth', 2, ...
        'color', 'black', 'displayname', 'Analytical');

    grid on;
    ylim([0, 20]);
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
    ylabel("$p_m(0, t)$", 'interpreter', 'latex');
    legend('location', 'northeast', 'interpreter', 'latex');

    %% Plot L2 norm error
%         axes('Position',[.34 .6 .3 .3])
%         box on
%         hold on;
%         ds_max = dns_ps(:, length(levels));
% 
% 
%         norms = zeros(length(levels) - 1, 1);
%         for level_idx = 1 : length(levels) - 1
%             level = levels(level_idx);
% 
%             % Determines L2-norm difference
%             diff = dns_ps(:, level_idx) - ds_max;
%             norms(level_idx) = sqrt(sum(diff.^2) / length(diff));
%         end
%         plot(levels(1 : end - 1), norms, 'color', 'black', 'linewidth', 2);
% 
%         sz = 100;
%         scatter(levels(1 : end - 1), norms, sz, colors, 'filled');
% 
%         set(gca, 'yscale', 'log');
%         xlim([min(levels) - 0.5, max(levels) - 0.5])
%         ylim([10^-3, 10^(-0.5)]);
%         xticks(levels(1 : end - 1));
%         yticks([10^-3, 10^-2, 10^-1]);
% 
%         grid on;
%         set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
%         xlabel("$m$", 'interpreter', 'latex');
%         ylabel("$||d_m - d_{14}||_2$", 'interpreter', 'latex');

    %% Overal figure settings
    x0=400;
    y0=400;
    height=800;
    width=600;

    set(gcf,'position',[x0,y0,width,height]);
    pause(1.5);

%     if imposed_coeff == "0"
%         figname = "dns_validation_figures/pressures/DNSPressure_Stationary";
%     else
%         figname = "dns_validation_figures/pressures/DNSPressure_Imposed";
%     end
%     set(gcf, 'Renderer', 'Painters');
%     exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
%     savefig(gcf, sprintf("%s.fig", figname));

end

    