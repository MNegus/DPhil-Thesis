%% FreeSurfaceComparison.m
% 
% 
% clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/FreeSurface/");
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
% Plate parameters
ALPHA = 2;
BETA = 0;
GAMMA = 20;
L = 2; 

% Computational parameters
DELTA_T = 1e-3;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
MAX_TIMESTEP = T_MAX / DELTA_T;

% DNS directories
dns_master_dir = "/home/michael/scratch/DPhil_DNS_Data";
stat_dir = sprintf("%s/Stationary_Plate/axi", dns_master_dir);
moving_dir = sprintf("%s/Moving_Plate/ALPHA-%g_BETA-%g_GAMMA-%g", dns_master_dir, ...
    ALPHA, BETA, GAMMA);

% Load DNS substrate solutions
moving_output_mat = readmatrix(sprintf("%s/cleaned_data/output.txt", moving_dir));
ts = moving_output_mat(:, 1) - IMPACT_TIME;
wsMoving = moving_output_mat(:, 6);
wsMovingFun = @(t) interp1(ts, wsMoving, t);

% Analytical parameters
epsilon = 1;

%% Load stationary analytical solutions
tMax = 1 / 3; % Max such that 1 = d(t)
tsStat = linspace(0, tMax, 1e3);

% Find substrate functions
StatSubstrateFunctions = platesubstratefunctions(tsStat, ...
    zeros(size(tsStat)), zeros(size(tsStat)), zeros(size(tsStat)), epsilon);

% Find ds
dsStat = StatSubstrateFunctions.d(tsStat);

%% Load moving analytical solutions
tMax = 1.5 * (1 / 3); % A bit after the stationary max value

% Solve plate equation
[tsComp, wsComp, w_tsComp, w_ttsComp] ...
    = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

% Find substrate functions
MovingSubstrateFunctions = platesubstratefunctions(tsComp, ...
    wsComp, w_tsComp, w_ttsComp, epsilon);

% Find where turnover point reaches 1
dsComp = MovingSubstrateFunctions.d(tsComp);
tIdxMaxComp = sum(dsComp <= 1);

% Restrict solutions temporally
tsComp = tsComp(1 : tIdxMaxComp);
wsComp = wsComp(1 : tIdxMaxComp);
w_tsComp = w_tsComp(1 : tIdxMaxComp);
w_ttsComp = w_ttsComp(1 : tIdxMaxComp);
dsComp = dsComp(1 : tIdxMaxComp);

%% Plot interfaces in time
% % Select timesteps to plot
% tsPlot = linspace(0.005, 0.2, 4);
% timesteps = floor((tsPlot + IMPACT_TIME) / DELTA_T);
% 
% % Set x limits
% xMins = [0.08, 0.44, 0.6, 0.7];
% plotWidths = [0.09, 0.09, 0.175, 0.29];
% yMins = -0.25 * plotWidths;
% 
% % Arrays for stationary and moving params
% dns_dirs = [stat_dir, moving_dir];
% function_structs = [StatSubstrateFunctions, MovingSubstrateFunctions];
% titleStrs = ["(a) Stationary substrate.", "(b) Moving substrate."];
% 
% % Loop over times
% for timestepIdx = 1 : length(timesteps)
% % for timestepIdx = 4
% 
%     timestep = timesteps(timestepIdx);
%     t = tsPlot(timestepIdx);
% 
%     tiledlayout(1, 2);
%     width = 3;
%     height = 3;
%     set(gcf,'units', 'inches', ...
%         'position',[0.5 * width, 0.5 * height, width, height]);
% 
%     % Loop over types
%     for typeIdx = 1 : 2
%         dns_dir = dns_dirs(typeIdx);
%         SubstrateFunctions = function_structs(typeIdx);
%         
%         nexttile;
%         hold on;
% 
%         %% Plot analytical solution
%         d = SubstrateFunctions.d(t);
%         J = SubstrateFunctions.J(t);
%         w = SubstrateFunctions.w(t);
%         xMaxUpper = 1;
%         xMaxLower = 1;
%     
%         % Load full composite solution
%         [xsTurnover, hsTurnover, ~, ~] ...
%             = outer_inner_jet_freesurface_composite(xMaxUpper, xMaxLower, ...
%             t, SubstrateFunctions);
% 
%         plot(xsTurnover, hsTurnover, 'linestyle', ':', ...
%             'color', redCol, 'linewidth', 1.25 * lineWidth, 'Displayname', 'Analytical');
%     
%         %% Plot DNS solution
%         interface_filename = sprintf("%s/interfaces/interface_%d.txt", ...
%             dns_dir, timestep);
%         transpose_coordinates = false;
%         
%         % Load interface points
%         [start_points, end_points] = ...
%             read_interface_points(interface_filename, transpose_coordinates);
%         
%         % Extract bulk droplet interface
%         tol = 1e-3;
%         [interface_start_points, interface_end_points] ...
%             = extract_interface(start_points, end_points, tol);
%         
%         % Determine substrate position
%         if typeIdx == 1
%             wVal = 0;
%         else
%             wVal = wsMovingFun(t);
%         end
% 
%         % Restrict vertical limits of points
%         yMax = yMins(timestepIdx) + plotWidths(timestepIdx);
%         keepIdxs = (interface_start_points(:, 1) < yMax + wVal) ...
%             & (interface_end_points(:, 1) < yMax + wVal);
%         interface_start_points = interface_start_points(keepIdxs, :);
%         interface_end_points = interface_end_points(keepIdxs, :);
% 
%         % Find line segments
%         xsStart = interface_start_points(:, 2);
%         ysStart = interface_start_points(:, 1) - wVal;
%     
%         xsEnd = interface_end_points(:, 2);
%         ysEnd = interface_end_points(:, 1) - wVal;     
% 
%         % Plot droplet interface
%         plot([xsStart'; xsEnd'], [ysStart'; ysEnd'], ...
%             'color', blueCol, 'linewidth', lineWidth);
% 
%         % Plot substrate
%         yline(-wVal, 'LineStyle', '--');
% 
%         % Figure properties
%         grid on;
%         box on;
%         xlabel("$r$");
%         ylabel("$z$");
%         xlim([xMins(timestepIdx), xMins(timestepIdx) + plotWidths(timestepIdx)]);
%         ylim([yMins(timestepIdx), yMins(timestepIdx) + plotWidths(timestepIdx)]);
%         pbaspect([1 1 1]);
% 
%     end
%     % Export figure
%     pause(0.5);
%     figname = append("PlateFigures/TurnoverInterfaceComparison_", num2str(t));
%     exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
%     % saveas(gcf, sprintf("%s.png", figname));
% end

%% Plot free surfaces in time
% Select timesteps to plot
tsPlot = linspace(0.005, 0.2, 4)
timesteps = floor((tsPlot + IMPACT_TIME) / DELTA_T);

% Variable colors
colorFreq = floor((length(cmap) / 3) / length(timesteps));

tiles = tiledlayout(2, 1);
width = 6;
height = 6;

% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);

% Set y limit to plot interface
yMax = 0.1;

% % Set x limits
% xWidth = 0.15;
% xStarts = [0.09, 0.43, 0.625, 0.72];

% Arrays for stationary and moving params
dns_dirs = [stat_dir, moving_dir];
function_structs = [StatSubstrateFunctions, MovingSubstrateFunctions];
titleStrs = ["(a) Stationary plate.", "(b) Moving plate."];


% Loop over types
for typeIdx = 1 : 2
    dns_dir = dns_dirs(typeIdx);
    SubstrateFunctions = function_structs(typeIdx);
    nexttile;
    hold on;

    for timestepIdx = 1 : length(timesteps)
        timestep = timesteps(timestepIdx);
        
        t = tsPlot(timestepIdx);

        % Set line colors
        AnalyticalLineColor = cmap(length(cmap) - (timestepIdx - 1) * colorFreq, :);
        DNSLineColor = cmap((timestepIdx - 1) * colorFreq + 1, :);
    
        %% Plot analytical solution
        d = SubstrateFunctions.d(t);
        J = SubstrateFunctions.J(t);
        w = SubstrateFunctions.w(t);
        xMaxUpper = 1;
        xMaxLower = 1;
    
        % Load full composite solution
        [xsTurnover, hsTurnover, ~, ~] ...
            = outer_inner_jet_freesurface_composite(xMaxUpper, xMaxLower, ...
            t, SubstrateFunctions);

        q = plot(xsTurnover, hsTurnover, 'linestyle', ':', ...
            'color', AnalyticalLineColor, 'linewidth', 1.25 * lineWidth, 'Displayname', 'Analytical');

        if timestepIdx == 1
            h(2) = q;
        end
%         scatter(d, (1 + 4 / pi) * J - w);
    
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
        
        % Determine substrate position
        if typeIdx == 1
            wVal = 0;
        else
            wVal = wsMovingFun(t);
%             wVal = 0;
        end

        % Restrict vertical limits of points
        keepIdxs = (interface_start_points(:, 1) < yMax + wVal) ...
            & (interface_end_points(:, 1) < yMax + wVal);
        interface_start_points = interface_start_points(keepIdxs, :);
        interface_end_points = interface_end_points(keepIdxs, :);

        % Find line segments
        xsStart = interface_start_points(:, 2);
        ysStart = interface_start_points(:, 1) - wVal;
    
        xsEnd = interface_end_points(:, 2);
        ysEnd = interface_end_points(:, 1) - wVal;     

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
        xlim([0, 1]);
        xlabel("$r$");
        ylabel("$z$");
        ylim([-0.05, yMax]);

    end

    %% Plot substrate
%     if typeIdx == 2
%         plot(dsComp, -wsComp, 'LineStyle', '--', 'Color', 0.75 * [1 1 1], ...
%             'LineWidth', lineWidth);
%     end

    %

    title(titleStrs(typeIdx), 'Fontsize', fontsize);
    set(gca, 'TitleFontSizeMultiplier', 1);

    drawnow;
    
end

%% Tiled plot properties
lh = legend([h(1), h(2)], ...
    ["DNS", "Analytical"], ...
    'NumColumns', 2);
% lh = legend('NumColumns', 2);
lh.Layout.Tile = 'North'; 

set(gcf, 'Renderer', 'Painters');
pause(0.5);

figname = "PlateFigures/TurnoverInterfaceComparison";
exportgraphics(tiles, sprintf("%s.png", figname), "Resolution", 300);
% saveas(gcf, sprintf("%s.png", figname));
% print(figname, '-dpng', '-r300');
