%% Free surface plot compare

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/FreeSurface");
addpath("../Analytical_Scripts/PlateSolution");
addpath("InterfaceAnalysis");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
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

%% DNS parameters
DELTA_T = 1e-3; % Time gap between interface files
IMPACT_TIME = 0.125; % Time of impact (theoretical)
DNS_T_MAX = 0.8; % Maximum DNS times
ts = - IMPACT_TIME : DELTA_T : DNS_T_MAX - IMPACT_TIME; % Time array
MAX_TIMESTEP = DNS_T_MAX / DELTA_T; % Index of maximum timestep
NO_TIMESTEPS = length(ts);

%% DNS directory name
master_dir = "/home/michael/scratch/DPhil_DNS_Data/Stationary_Plate/axi";

%% Load analytical solutions (stationary case)
epsilon = 0.1;
tsAnalytical = (0 : DELTA_T : DNS_T_MAX - IMPACT_TIME) / epsilon^2;
ws = zeros(size(tsAnalytical));
w_ts = zeros(size(tsAnalytical));
w_tts = zeros(size(tsAnalytical));
SubstrateFunctions = platesubstratefunctions(tsAnalytical, ...
    ws, w_ts, w_tts, epsilon);

%% Plot interface in time
figure(1);
% for k = 1 : NO_TIMESTEPS
for k = 150 : 200
    t = ts(k); % Time variable
    
    %% Plot DNS solution
    interface_filename = sprintf("%s/interfaces/interface_%d.txt", ...
        master_dir, k);
    transpose_coordinates = false;
    
    % Load interface points
    [start_points, end_points] = ...
        read_interface_points(interface_filename, transpose_coordinates);
    
    % Extract bulk droplet interface
    tol = 1e-3;
    [interface_start_points, interface_end_points] ...
        = extract_interface(start_points, end_points, tol);
    
%     % Find line segment centres
%     seg_centres = 0.5 * (interface_start_points + interface_end_points);
%     [hsDNS, sort_idxs] = sort(seg_centres(:, 1));
%     xsDNS = seg_centres(sort_idxs, 2);
    allPoints = [interface_start_points, interface_end_points]
    xsDNS = allPoints(:, 2);
    hsDNS = allPoints(:, 1);
    
    % Plot droplet interface
%     plot(xsDNS, hsDNS);
    scatter(xsDNS, hsDNS);
    
    %% Plot analytical solution
    if (t > 0)
        tAnalytical = t / epsilon^2;
        d = SubstrateFunctions.d(tAnalytical);
        xMaxUpper = 1;
        xMaxLower = 1.5 * epsilon * d;
        % Load full composite solution
        [xsTurnover, hsTurnover, ~, ~] ...
            = outer_inner_jet_freesurface_composite(xMaxUpper, xMaxLower, ...
            tAnalytical, SubstrateFunctions);
        hold on;
        plot(xsTurnover, hsTurnover);
%         plot(xsFull, hsFull);
        hold off;
    end
    
    %% Figure options
    xlim([epsilon * d - 0.1, epsilon * d + 0.1]);
    ylim([0, 0.1]);
%     xlim([0, 2]);
%     ylim([0, 2]);
    pbaspect([1 1 1]);
    pause(0.1);
end