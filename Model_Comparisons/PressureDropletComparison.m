%% Whole droplet pressure

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Forces/");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;
colormap(cmap);

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
DELTA_T = 1e-2;
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
stat_dir = sprintf("%s/Stationary_Plate/axi", dns_master_dir);
moving_dir = sprintf("%s/Moving_Plate/ALPHA-2_BETA-0_GAMMA-20", dns_master_dir);

%% Load DNS substrate position


%% Plot pressure at different timesteps
timesteps = 80;
pLimits = [0.9];

% Arrays for stationary and moving params
dns_dirs = [stat_dir, moving_dir];
% function_structs = [StatSubstrateFunctions, MovingSubstrateFunctions];
titleStrs = ["(a) Stationary substrate.", "(b) Moving substrate."];


figure(1);
hold on;
for timestep = timesteps
    
    % Loop over types
    for typeIdx = 1 : 2
        % Load in the field struct
        loadObj = load(sprintf("%s/coarsenedOutputs/fieldStruct_%d.mat", ...
            dns_dirs(typeIdx), timestep));
        fieldStruct = loadObj.fieldStruct;

        % 
        surf((2 * (typeIdx - 1) - 1) * fieldStruct.Y, fieldStruct.X, fieldStruct.ps, 'EdgeColor','none');
            

    end
    
    view(2);
    pbaspect([2, 1, 1]);
end

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
