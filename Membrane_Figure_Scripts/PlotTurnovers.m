close all;
clear;

addpath("FiniteDifference/");
addpath("FiniteDifference/PressuresFD/");
addpath("NormalModes/");

parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";


%% Parameters
EPSILON = 1;
L = 16;
T_MAX = 0.25;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 10924;

xs = linspace(0, L, N_MEMBRANE)';

ts_analytical = 0 : DELTA_T : T_MAX;
size(ts_analytical)

%% Loop over GAMMAS
BETA = 0;
ALPHA = 1.002;
GAMMAS = 133.6 * 10.^linspace(-2, 2, 101);
gamma_vary_dir = sprintf("%s/GAMMA_varying", parent_dir);

figure(1);
hold on;

figure(2);
hold on;

for GAMMA = GAMMAS(1 : end)
    %% GAMMA directory
    param_dir = sprintf("%s/GAMMA_%g", gamma_vary_dir, GAMMA)
    
    %% Finite differences
%     fd_data_dir = sprintf("%s/FiniteDifference", param_dir)
%     mkdir(fd_data_dir);
%     mkdir(append(fd_data_dir, "/composite"));
%     SaveFDSolution(fd_data_dir, ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, "composite");

    %% Norma modes
    nm_data_dir = sprintf("%s/NormalModes", param_dir);
    
    % Load struct
    SolStruct = load(append(nm_data_dir, "/SolStruct.mat")).SolStruct;

    figure(1);
    plot(SolStruct.ts, SolStruct.ds);
    drawnow;

    %% Membrane solution
    k = 2500;
    [ws, w_ts, ps] ...
        = MembraneSolutionNM(xs, SolStruct.as(k, :), ...
        SolStruct.a_ts(k, :), SolStruct.q_ts(k, :), ...
        SolStruct.ds(k), L, SolStruct.N, EPSILON);
    figure(2);
    plot(xs, -ws);
end

figure(1);
plot(ts_analytical, 2 * sqrt(ts_analytical), 'color', 'black', 'LineStyle','--')