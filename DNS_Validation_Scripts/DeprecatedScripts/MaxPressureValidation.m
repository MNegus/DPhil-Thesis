%% MaxPressureValidation.m
% Validates the pointwise maximum pressure along the substrate

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/Imposed/");
addpath("../Analytical_Scripts/Imposed/Pressures");

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

fontsize = 22;

%% Load in analytical parameters
analyticalType = "curvedDNS";
[epsilon, q, omega, p, L] = substrateparameters(analyticalType);

%% Computational parameters
DELTA_T = 1e-4; % Timestep
IMPACT_TIME = 0.125; % Time of impact
% T_MAX = 0.25 + IMPACT_TIME; % Maximum time to plot
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME; % Time array
% tsAnalytical = DELTA_T : DELTA_T : T_MAX - IMPACT_TIME; % Analytical time array
tsAnalytical = 1e-9 : DELTA_T : 0.25; % Analytical time array
MAX_TIMESTEP = T_MAX / DELTA_T; % Maximum timestep
NO_TIMESTEPS = length(ts); % Number of timesteps
freq = 100; % Frequency of scatter

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
type = ["imposed_0.05_omega_4"];
levels = [10, 11, 12, 13, 14];

%% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels - 1, 3);

for q = 1 : no_levels - 1
    colors(q, :) = cmap(color_idxs(q), :);
end


%% Plot all cases
imposed_coeffs = ["0"; "0.05"]; % Coefficients of imposed solutions
for imposedIdx = 1 : 2
    
    imposed_coeff = imposed_coeffs(imposedIdx); % Save coefficient
    
    % Set options
    if imposed_coeff == "0"
        analyticalType = "stationary";
%         parent_dir = "/scratch/negus/DNS_Chapter_validation/imposed_coeff_0";
    else
        analyticalType = "curvedDNS";
    end
    
    parent_dir = sprintf("%s/%s_maxlevel_validation/imposed_coeff_%s", ...
            master_dir, type, imposed_coeff);

    %% Analytical solutions
    % Load in substrate functions
    [epsilon, ~, ~, ~, ~] = substrateparameters(analyticalType)
    SubstrateFunctions = substratefunctions(analyticalType, "2D");
    
    % Determine maximum pressure
    [pMaxsAnalytical, ~] = pressuremax(tsAnalytical / epsilon^2, SubstrateFunctions);
    
    % Determine pressure at the origin
    p0sAnalytical = zeros(size(tsAnalytical));
    for k = 1 : length(tsAnalytical)
        t = tsAnalytical(k);
        p0sAnalytical(k) = outerpressure(0, t / epsilon^2, SubstrateFunctions);
    end
    
    %% Create figures
    % Pressure max figure
%     maxIdx = imposedIdx;
%     originIdx = imposedIdx + length(imposed_coeffs);
    
%     close(figure(maxIdx));
%     figure(maxIdx);
%     hold on;
    
    originIdx = imposedIdx;
%     close(figure(originIdx));
    figure(originIdx);
    hold on;

    % Matrices to hold DNS pressure values
    pMaxsDNS = zeros(NO_TIMESTEPS, length(levels));
    p0sDNS = zeros(NO_TIMESTEPS, length(levels));
    
    %% Loop over the levels
    for level_idx = 1 : length(levels)
        
        % Find computational directory
        level = levels(level_idx);
        level_dir = sprintf("%s/max_level_%d", parent_dir, level);
        
        % Find origin pressure from output.txt
        output_mat = readmatrix(sprintf("%s/output.txt", level_dir));
        tsDNS = output_mat(:, 1);
        p0sDNS(:, level_idx) = output_mat(1 : NO_TIMESTEPS, 3);
        
%         % Loop over timesteps
%         for timestep = 1 : freq : NO_TIMESTEPS
%             timestep
%             
%             % Load in boundary value matrix
%             membrane_mat = readmatrix(sprintf("%s/boundary_outputs/boundary_output_%d.txt", level_dir, timestep - 1));
%             
%             % Sort matrix
%             xs = membrane_mat(:, 1);
%             [sort_xs, sort_idxs] = sort(xs);
%             ps = membrane_mat(sort_idxs, 2);
%             
%             % Find maximum pressure
%             pMaxsDNS(timestep, level_idx) = max(ps);
%             
%             % Find origin pressure
%             p0sDNS(timestep, level_idx) = ps(1);
%         end

        % Plot maximum pressure
%         figure(maxIdx);
%         scatter(ts(1 : freq : end), pMaxsDNS(1 : freq : end, level_idx), [], cmap(color_idxs(level_idx), :), ...
%             'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
%         drawnow;
        
        % Plot origin pressure
        figure(originIdx);
        scatter(ts(1 : freq : end), p0sDNS(1 : freq : end, level_idx), [], cmap(color_idxs(level_idx), :), ...
            'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
%         scatter(ts(1 : freq : end), p0sDNS(1 : freq : end, level_idx), [], ...
%             'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
        drawnow;
    end
    
    %% Plot analytical solutions
%     figure(maxIdx);
%     plot(tsAnalytical, pMaxsAnalytical, 'linestyle', '--', 'linewidth', 2, ...
%         'color', 'black', 'displayname', 'Analytical');
%     ylabel("max($p_m(x, t)$)", 'interpreter', 'latex');
    
    figure(originIdx);
    plot(tsAnalytical, p0sAnalytical, 'linestyle', '--', 'linewidth', 2, ...
        'color', 'black', 'displayname', 'Analytical');

    ylabel("$p_m(0, t)$", 'interpreter', 'latex');

    % Figure settings
%     for idx = [maxIdx, originIdx]
    for idx = originIdx
        figure(idx);
        grid on;
        ylim([0, 20]);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
        xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
        legend('location', 'northeast', 'interpreter', 'latex');
        
        x0=400;
        y0=400;
        height=800;
        width=600;

        set(gcf,'position',[x0,y0,width,height]);
        pause(1.5);
    end
    

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


%     if imposed_coeff == "0"
%         figname = "dns_validation_figures/pressures/DNSPressure_Stationary";
%     else
%         figname = "dns_validation_figures/pressures/DNSPressure_Imposed";
%     end
%     set(gcf, 'Renderer', 'Painters');
%     exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
%     savefig(gcf, sprintf("%s.fig", figname));

end