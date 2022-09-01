%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD

addpath("PressuresFD");
addpath("../");

parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";


%% Parameters
EPSILON = 1;
ALPHA = 0.0017; BETA = 0.681; GAMMA = 0.003; 
L = 2;
T_MAX = 0.575;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 10924;

%% Finite differences
fd_data_dir = sprintf("%s/FiniteDifference", parent_dir);
SaveFDSolution(fd_data_dir, ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, "composite");