%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD

addpath("FiniteDifference/");
addpath("FiniteDifference/PressuresFD/");
addpath("NormalModes/");

parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";

%% Parameters
EPSILON = 1;
L = 16;
T_MAX = 0.35;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 10924;

xs = linspace(0, L, N_MEMBRANE)';

ts_analytical = 0 : DELTA_T : T_MAX;

%% Rubber sheet case
ALPHA = 1.1; BETA = 0; GAMMA = 668.0;
rubber_sheet_dir = sprintf("%s/RubberSheet", parent_dir);
mkdir(rubber_sheet_dir);

% Finite difference
fd_data_dir = sprintf("%s/FiniteDifference", rubber_sheet_dir);
mkdir(fd_data_dir);
mkdir(append(fd_data_dir, "/composite"));
SaveFDSolution(fd_data_dir, ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, "composite");

% Normal modes
nm_data_dir = sprintf("%s/NormalModes", rubber_sheet_dir);
mkdir(nm_data_dir);
SaveValidatedNMSolution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX, DELTA_T);



%% Loop over GAMMAS
BETA = 0;
ALPHA = 1.002;
GAMMAS = 133.6 * 10.^linspace(-2, 2, 101);
gamma_vary_dir = sprintf("%s/GAMMA_varying", parent_dir);
mkdir(gamma_vary_dir)
for GAMMA = GAMMAS(1 : end)
    %% GAMMA directory
    param_dir = sprintf("%s/GAMMA_%g", gamma_vary_dir, GAMMA)
    mkdir(param_dir)
    
    %% Finite differences
%     fd_data_dir = sprintf("%s/FiniteDifference", param_dir)
%     mkdir(fd_data_dir);
%     mkdir(append(fd_data_dir, "/composite"));
%     SaveFDSolution(fd_data_dir, ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, "composite");

    %% Normal modes
    nm_data_dir = sprintf("%s/NormalModes", param_dir);
    mkdir(nm_data_dir);
    SaveValidatedNMSolution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX, DELTA_T);

end

% %% Finite differences
% fd_data_dir = sprintf("%s/FiniteDifference", parent_dir);
% SaveFDSolution(fd_data_dir, ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, "composite");
% 
% 
% %% Validated NM
% nm_data_dir = sprintf("%s/NormalModes", parent_dir);
% SaveValidatedNMSolution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX, DELTA_T);