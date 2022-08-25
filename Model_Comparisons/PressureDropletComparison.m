%% Whole droplet pressure

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Forces/");
addpath("../DNS_Post-Processing/InterfaceAnalysis/");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;

redCol = cmap(end, :);
blueCol = cmap(1, :);

%% Figure options
fontsize = 8.5;
lineWidth = 1;
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize', fontsize);
set(0, 'DefaultTextFontSize', fontsize);
set(0,'defaultLegendFontSize', fontsize, 'DefaultLegendFontSizeMode','manual');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(0, 'DefaultLegendInterpreter', 'latex');

%% Parameters
% Computational parameters
DELTA_T = 1e-2;
INTERFACE_DELTA_T = 1e-3;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
MAX_TIMESTEP = T_MAX / DELTA_T;

% Master DNS dir
dns_master_dir = "/home/michael/scratch/DPhil_DNS_Data/";

% Analytical parameters
epsilon = 1;
tMax = 1.5 * (1 / 3); % A bit after the stationary max value

% Plate parameters
ALPHA = 2;
BETA = 0;
GAMMA = 20;

%% DNS directories
stat_dir = sprintf("%s/Stationary_Plate_Output/axi", dns_master_dir);
moving_dir = sprintf("%s/Moving_Plate_Output/", dns_master_dir);

%% Load DNS substrate position
moving_output_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", moving_dir));
ts = moving_output_mat(:, 1) - IMPACT_TIME;
wsMoving = moving_output_mat(:, 6);
wsMovingFun = @(t) interp1(ts, wsMoving, t);

%% Plot pressure at different timesteps
timesteps = [13, 20, 50, 80];
pLimits = [13, 4.9, 1.25, 0.9];

% Arrays for stationary and moving params
dns_dirs = [stat_dir, moving_dir];
% function_structs = [StatSubstrateFunctions, MovingSubstrateFunctions];
titleStrs = ["(a) Stationary substrate.", "(b) Moving substrate."];

% tiledlayout(2, 2);
% colormap(cmap);
% width = 4;
% height = 5;
% 
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);

for timestepIdx = 1 : length(timesteps)
%     nexttile;

    figure(timestepIdx);
    colormap(cmap);
    width = 2.5;
    height = 2;
    
    set(gcf,'units', 'inches', ...
        'position',[0.5 * width, 0.5 * height, width, height]);

    
    hold on;

    % Timestep 
    timestep = timesteps(timestepIdx);

    % Save time
    t = DELTA_T * timestep - IMPACT_TIME

    % Loop over types
    for typeIdx = 1 : 2
        % DNS directory
        dns_dir = dns_dirs(typeIdx);

        %% Load in the field struct
        loadObj = load(sprintf("%s/coarsenedOutputs/fieldStruct_%d.mat", ...
            dns_dir, timestep));
        fieldStruct = loadObj.fieldStruct;

        %% Load substrate solution
        if typeIdx == 1
            wVal = 0;
        else
            wVal = wsMovingFun(t);
        end

        % Load in mesh grid (need to transpose data)
        Xs = fieldStruct.Y;
        Ys = fieldStruct.X - wVal;
        ps = fieldStruct.ps;

        % Factor in front of Xs to reflect the plot
        refFactor = (2 * (typeIdx - 1) - 1);

        %% Surf plot of pressure
%         surf(refFactor * Xs, Ys, ps, 'EdgeColor','none');
        s = pcolor(refFactor * Xs, Ys, ps);
        s.EdgeColor = 'none';
        clim([0, pLimits(timestepIdx)]);
%         view(2);
        
        %% Interface plot
        interface_timestep = floor(timestep * DELTA_T / INTERFACE_DELTA_T);
        
        interface_filename = sprintf("%s/interfaces/interface_%d.txt", ...
            dns_dir, interface_timestep);
        transpose_coordinates = false;
        
        % Load interface points
        [start_points, end_points] = ...
            read_interface_points(interface_filename, transpose_coordinates);
        
        % Extract bulk droplet interface
        tol = 1e-3;
        [interface_start_points, interface_end_points] ...
            = extract_interface(start_points, end_points, tol);

        % Restrict vertical limits of points
        yMax = 2;
        keepIdxs = (interface_start_points(:, 1) < yMax + wVal) ...
            & (interface_end_points(:, 1) < yMax + wVal);
        interface_start_points = interface_start_points(keepIdxs, :);
        interface_end_points = interface_end_points(keepIdxs, :);

        % Find line segments
        freq = 1;
        xsStart = interface_start_points(1 : freq : end, 2);
        ysStart = interface_start_points(1 : freq : end, 1) - wVal;
    
        xsEnd = interface_end_points(1 : freq : end, 2);
        ysEnd = interface_end_points(1 : freq : end, 1) - wVal;     

        % Plot droplet interface
%         zPos = pLimits(timestepIdx) * ones(size([ysStart'; ysEnd']));
        plot([refFactor * xsStart'; refFactor * xsEnd'], [ysStart'; ysEnd'], ...
             'color', 'black', 'linewidth', lineWidth);

        %% Plot substrate position
        xs = [0, 2];
        ws = wVal * ones(size(xs));
        plot(refFactor * xs, -ws, 'Linestyle', '--', ...
            'Color', 0.75 * [1 1 1], 'LineWidth', lineWidth);
    end

    %% Add white dividing line
    xline(0, 'color', 'white');

    %% Add color bar
    cb = colorbar('Location', 'Northoutside');
    cb.Label.Interpreter = 'latex';
    cb.TickLabelInterpreter = 'latex';
    cb.Label.String = '$p$';

    %% Set figure options
    xMin = -2;
    xMax = 2;
    xlim([xMin, xMax]);
    xlabel("$r$");
    
    yMin = -0.25;
    yMax = 2.1;
    ylim([yMin, yMax]);
    ylabel("$z$");

    pbaspect([xMax - xMin, yMax - yMin, 1]);

    pause(0.5);
    
    figname = "PlateFigures/DNS_Pressure";
    exportgraphics(gcf, sprintf("%s_%g.png", figname, t), "Resolution", 300);
end

% figname = "PlateFigures/DNS_Pressure";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%%
% data = readmatrix("field_output_13.txt");
% 
% 
% %%
% N = 2731;
% 
% psFull = data(:, 3);
% 
% psMesh = reshape(psFull, [N, N]);
% 
% %% 
% xsFull = linspace(0, 2, N);
% ysFull = linspace(0, 2, N);
% [XFull, YFull] = meshgrid(xsFull, ysFull);
% 
% NCoarse = 256;
% xs = linspace(0, 2, NCoarse);
% ys = linspace(0, 1, NCoarse);
% [X, Y] = meshgrid(xs, ys);
% psCoarse = interp2(XFull, YFull, psMesh, X, Y);
% 
% %% contourf
% shading interp
% surf(X, Y, psCoarse, 'EdgeColor','none');
