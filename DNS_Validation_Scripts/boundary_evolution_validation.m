%% pressure_evolution_validation.m
% Validates the pressure evolution in the DNS by comparing its value 
% for multiple levels

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/Imposed/");
addpath("../Analytical_Scripts/Imposed/Pressures");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;

fontsize = 22;

%% Parameters
L = 2;
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
IMPACT_TIMESTEP = IMPACT_TIME / DELTA_T;
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);
freq = 5; % Frequency

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
type = ["imposed_0.05_omega_4"];
levels = 11:12;

%% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels - 1, 3);

for q = 1 : no_levels - 1
    colors(q, :) = cmap(color_idxs(q), :);
end

%% Plot all cases

imposed_coeffs = ["0", "0.05"];

for imposedIdx = 2
    imposed_coeff = imposed_coeffs(imposedIdx);
   
    
    %% Set options
    if imposed_coeff == "0"
        analyticalType = "stationary";
%         parent_dir = "/scratch/negus/DNS_Chapter_validation/imposed_coeff_0";
    else
        analyticalType = "curvedDNS";
%         analyticalType = "stationary";
%         parent_dir = sprintf("%s/%s_maxlevel_validation/imposed_coeff_%s", ...
%             master_dir, type, imposed_coeff);
    end
    
    parent_dir = sprintf("%s/%s_maxlevel_validation/imposed_coeff_%s", ...
            master_dir, type, imposed_coeff);
    %% Analytical solutions
    % Load in substrate functions
    [epsilon, ~, ~, ~, ~] = substrateparameters(analyticalType);
    SubstrateFunctions = substratefunctions(analyticalType, "2D");
    
    % Analytical parameters
    xsAnalytical = linspace(0, L, 1e3);
    
    %% Video setup
    vidName = sprintf("dns_validation_figures/imposed_coeff_%s.mp4", imposed_coeffs(imposedIdx));
    writerObj = VideoWriter(vidName);
    writerObj.FrameRate = 10;
    open(writerObj);
    frameCounter = 0;

    
    layout = tiledlayout(2, 1);
    
    %% Loop over time
    for timestep = 1300 : 10: 1700
        tAnalytical = DELTA_T * timestep - IMPACT_TIME;

        %% Loop over the levels
        for level_idx = 1 : length(levels)
            
            % Level directory
            level = levels(level_idx);
            level_dir = sprintf("%s/max_level_%d", parent_dir, level);
            
            % Load DNS solutions
            membrane_mat = readmatrix(sprintf("%s/boundary_outputs/boundary_output_%d.txt", level_dir, timestep));
            xs = membrane_mat(:, 1);
            [sort_xs, sort_idxs] = sort(xs);
            ps = membrane_mat(sort_idxs, 2);
            vs = membrane_mat(sort_idxs, 4);
            
            nexttile(1);
            plot(sort_xs, vs, 'color', cmap(color_idxs(level_idx), :), ...
                'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
            hold on;
            
            
            nexttile(2);
            plot(sort_xs, ps, 'color', cmap(color_idxs(level_idx), :), ...
                'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
            hold on;
            
        
            
        end
        
        %% Plot analytical solutions
        % Substrate solution
        nexttile(1);
        w_tsAnalytical = SubstrateFunctions.w_t(xsAnalytical, tAnalytical / epsilon^2);
        plot(xsAnalytical, -w_tsAnalytical, 'linewidth', 2, ...
        'color', 'black', 'displayname', 'Analytical');
        
        % Pressure
        nexttile(2);
        [psAnalytical, ~, ~] ...
            = substratepressure(xsAnalytical, tAnalytical / epsilon^2, SubstrateFunctions);
        plot(xsAnalytical, psAnalytical, 'linewidth', 2, ...
        'color', 'black', 'displayname', 'Analytical');
    
        [pMaxAnalytical, xMaxAnalytical] = pressuremax(tAnalytical / epsilon^2, SubstrateFunctions);
        scatter(xMaxAnalytical, pMaxAnalytical, 'filled');
    
%         ylim([-0.1 * max(psAnalytical), 1.5 * max(psAnalytical)]);
        ylim([-0.1, 20]);

        %% Velocity figure settings
        nexttile(1);
        hold off;
        grid on;
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
        xlabel("$x$", 'interpreter', 'latex', "Fontsize", fontsize);
        ylabel("$v(x, 0, t)$", 'interpreter', 'latex');
        legend('location', 'northeast', 'interpreter', 'latex');
        xlim([0, 2]);
        ylim([-0.5, 0.5]); 
        
        %% Pressure figure settings
        nexttile(2);
        hold off;
        grid on;
        xlim([0, 2]);
        set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
        xlabel("$x$", 'interpreter', 'latex', "Fontsize", fontsize);
        ylabel("$p(x, 0, t)$", 'interpreter', 'latex');
        legend('location', 'northeast', 'interpreter', 'latex');
        
        
        
        %% Overal figure settings
        
        title(layout,"$t = $" + num2str(tAnalytical), 'interpreter', 'latex', ...
            'Fontsize', 24);
        x0=400;
        y0=400;
        height=800;
        width=800;

        set(gcf,'position',[x0,y0,width,height]);
        drawnow;
        
        if frameCounter > 5
            frame = getframe(gcf);
            writeVideo(writerObj, frame);
        end
        
        frameCounter = frameCounter + 1;

        
%         pause(1.5);

    end
    
    close(writerObj);
end
    