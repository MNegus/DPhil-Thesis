%% FreeSurfaceComparison.m
% 
% 
% clear;
close all;

% Adds analytical scripts to path
addpath("../DNS_Post-Processing/InterfaceAnalysis/");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;

redCol = cmap(end, :);
blueCol = cmap(1, :);

%% Figure options
fontsize = 10;
lineWidth = 1.25;
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize', fontsize);
set(0, 'DefaultTextFontSize', fontsize);
set(0,'defaultLegendFontSize', fontsize, 'DefaultLegendFontSizeMode','manual');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(0, 'DefaultLegendInterpreter', 'latex');

%% Parameters

% Computational parameters
DELTA_T = 1e-3;
IMPACT_TIME = 0.125;
T_MAX = 0.7;
MAX_TIMESTEP = T_MAX / DELTA_T;

% DNS directories
dns_master_dir = "/home/michael/scratch/DPhil_DNS_Data";
stat_dir = sprintf("%s/Stationary_Membrane", dns_master_dir);
moving_dir = sprintf("%s/RubberRuns/ALPHA_1.1_GAMMA_668", dns_master_dir);

% Analytical parameters
epsilon = 1;


%% Plot free surfaces in time (t up to 0.2)
% Select timesteps to plot
tsPlot = linspace(0.005, 0.575, 5)
timesteps = floor((tsPlot + IMPACT_TIME) / DELTA_T);

% Variable colors
colorFreq = floor((length(cmap) / 3) / length(timesteps));

tiles = tiledlayout(2, 1);
width = 6;
height = 6;

% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);

% Set y limit to plot interface
yMax = 0.5;

% % Set x limits
% xWidth = 0.15;
% xStarts = [0.09, 0.43, 0.625, 0.72];

% Arrays for stationary and moving params
dns_dirs = [stat_dir, moving_dir];
titleStrs = ["(a) Stationary membrane.", "(b) Moving membrane."];


% Loop over types
for typeIdx = 1 : 2
    dns_dir = dns_dirs(typeIdx);
    nexttile;
    hold on;

    for timestepIdx = 1 : length(timesteps)
        timestep = timesteps(timestepIdx);
        
        t = tsPlot(timestepIdx);

        % Set line colors
        DNSLineColor = cmap((timestepIdx - 1) * colorFreq + 1, :);
    
        %% Plot DNS solution
        interface_filename = sprintf("%s/interfaces/interface_%d.txt", ...
            dns_dir, timestep);
        transpose_coordinates = false;
        
        % Load interface points
        [start_points, end_points] = ...
            read_interface_points(interface_filename, transpose_coordinates);
        
        % Extract bulk droplet interface
        tol = 1e-3;
        [interface_start_points, interface_end_points] ...
            = extract_interface(start_points, end_points, tol);
        

        % Restrict vertical limits of points
        keepIdxs = (interface_start_points(:, 1) < yMax) ...
            & (interface_end_points(:, 1) < yMax);
        interface_start_points = interface_start_points(keepIdxs, :);
        interface_end_points = interface_end_points(keepIdxs, :);

        % Find line segments
        xsStart = interface_start_points(:, 2);
        ysStart = interface_start_points(:, 1);
    
        xsEnd = interface_end_points(:, 2);
        ysEnd = interface_end_points(:, 1);     

        % Plot droplet interface
        plot([xsStart'; xsEnd'], [ysStart'; ysEnd'], ...
            'color', DNSLineColor, 'linewidth', lineWidth);

        % Plot first element for legend
        if timestepIdx == 1
            h(1) = plot([xsStart(1)'; xsEnd(1)'], [ysStart(1)'; ysEnd(1)'], ...
            'color', DNSLineColor, 'linewidth', lineWidth);
        end

        % Figure properties
        grid on;
        box on;
        xlim([0, 2.5]);
        xlabel("$x$");
        ylabel("$z$");
        ylim([-0.05, yMax]);

    end
    
    title(titleStrs(typeIdx), 'Fontsize', fontsize);
    set(gca, 'TitleFontSizeMultiplier', 1);

    drawnow;
    
end

set(gcf, 'Renderer', 'Painters');
pause(0.5);

figname = "MembraneFigures/MembraneInterfaceComparison";
exportgraphics(tiles, sprintf("%s.png", figname), "Resolution", 300);
