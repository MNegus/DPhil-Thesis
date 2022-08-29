%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD

addpath("PressuresFD");
addpath("../");

parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";


%% Parameters
EPSILON = 1;
ALPHA = 1.002; BETA = 0.0; GAMMA = 133.6; 
L = 16;
T_MAX = 0.575;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 10924;

%% Validated NM
[N, delta_d, ts, ds, as, a_ts, a_tts, q_ts] ...
    = ValidatedNMSolution(ALPHA, BETA, GAMMA, EPSILON, L, T_MAX, DELTA_T);

%% 