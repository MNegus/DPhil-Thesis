%% OriginPressureValidation.m
% Validates the force in the DNS by comparing it for multiple levels

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/Imposed/");
addpath("../Analytical_Scripts/Imposed/Pressures");

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
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
L = 2;
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);
freq = 20; % Frequency of scatter

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
type = ["imposed_0.05_omega_4"];
level = 13;
imposed_coeffs = ["0", "0.05"]; % Coefficients of imposed substrate

%% Analytical solutions
xsAnalytical = linspace(0, L, 1e3); % Array of x values

% Stationary substrate solutions
[epsilonStationary, ~, ~, ~, ~] = substrateparameters("stationary");
StationaryFunctions = substratefunctions("stationary", "2D");

% Curved substrate solutions
[epsilonMoving, ~, ~, ~, ~] = substrateparameters("curvedDNS");
MovingFunctions = substratefunctions("curvedDNS", "2D");

%% Plot all cases
% Select timesteps to plot
timesteps = [2000, 3000];

% Set up tiled figure
tileFig = tiledlayout(1, 2, 'Tilespacing', 'compact', 'padding', 'compact');

% Loop over time
for timestepIdx = 1 : 2
    timestep = timesteps(timestepIdx);
    
    tAnalytical = DELTA_T * timestep - IMPACT_TIME; % Analytical time
    
    % Move to relevant tile
    nexttile(timestepIdx);
    hold on;
    
    for imposedIdx = 1 : 2
        imposed_coeff = imposed_coeffs(imposedIdx); % Save coefficient
        
        % Parent directory
        parent_dir = sprintf("%s/%s_maxlevel_validation/imposed_coeff_%s", ...
            master_dir, type, imposed_coeff);
        
        % Level directory
        level_dir = sprintf("%s/max_level_%d", parent_dir, level); 
        
        % Load DNS solutions
        membrane_mat = readmatrix(sprintf("%s/boundary_outputs/boundary_output_%d.txt", level_dir, timestep));
        xs = membrane_mat(:, 1);
        [sort_xs, sort_idxs] = sort(xs);
        ps = membrane_mat(sort_idxs, 2);
        
        % Save displayname
        if imposedIdx == 1
            displayName = "DNS: Stationary substrate";
            colorIdx = floor(0.1 * length(cmap));
        else
            displayName = "DNS: Moving substrate";
            colorIdx = ceil(0.9 * length(cmap));
        end
        
        % Plot pressure
        h(imposedIdx) = plot(sort_xs, ps, 'color', cmap(colorIdx, :), ...
            'linewidth', lineWidth, 'Displayname', displayName);
        hold on;
    end
    
    % Plot stationary analytical solutions
    [psStationary, ~, ~] ...
            = substratepressure(xsAnalytical, tAnalytical / epsilonStationary^2, ...
            StationaryFunctions);
    h(3) = plot(xsAnalytical, psStationary, 'linestyle', '--', ...
        'color', 'black', 'linewidth', lineWidth, ...
        'Displayname', "Analytical: Stationary substrate");
    
    % Plot moving analytical solutions
    [psMoving, ~, ~] ...
            = substratepressure(xsAnalytical, tAnalytical / epsilonMoving^2, ...
            MovingFunctions);
    h(4) = plot(xsAnalytical, psMoving, 'linestyle', ':', ...
        'color', 'black', 'linewidth', 1.5 * lineWidth, ...
        'Displayname', "Analytical: Moving substrate");
    
    hold off;
    
    % Tile figure options
    box on;
    grid on;
    xlim([0, 1.25]);
    xticks([0, 0.25, 0.5, 0.75, 1, 1.25]);
    
    if timestepIdx == 1
        ylim([-1, 7]);
    else
        ylim([-0.5, 3.5]);
    end
    
    xlabel("$x$");
    ylabel("$p(x, t)$");
    
    titleStr = "$t = " + num2str(tAnalytical) + "$";
    title(titleStr);
    
    
end

%% Whole figure options
lh = legend(h(1:4), 'Numcolumns', 2);
lh.Layout.Tile = 'South'; 

% Set size in inches
width = 6;
height = 3.65;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

%% Save figure
set(gcf, 'Renderer', 'Painters');
pause(1.5);

% Figure name
figname = "dns_validation_figures/pressures/DNSPressureEvolution";

% Export figure
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
savefig(gcf, sprintf("%s.fig", figname));
