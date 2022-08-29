%% MembraneDisplacementComparison
% 

clear;
close all;

% Add paths
addpath("FiniteDifference/");
addpath("FiniteDifference/PressuresFD/");
addpath("NormalModes/");

% Parent directory where data is stored
parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";

param_dir = append(parent_dir, "RubberSheet");

%% Parameters
EPSILON = 1;
L = 16;
T_MAX = 0.35;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 10924;

xs = linspace(0, L, N_MEMBRANE)';

ts_analytical = 0 : DELTA_T : T_MAX;

%% Load normal modes solutions
SolStruct = load(sprintf("%s/NormalModes/SolStruct.mat", param_dir)).SolStruct;
as = SolStruct.as;
a_ts = SolStruct.a_ts;
q_ts = SolStruct.q_ts;
ds = SolStruct.ds;
N_M = SolStruct.N;

%% Plot membrane in time
xMax = L;
tiledlayout(3, 1);

for k = 2 : 10 : length(ts_analytical)
    t = ts_analytical(k)

    %% Normal modes solution
    [ws, w_ts, ps] ...
        = MembraneSolutionNM(xs, as(k, :), ...
        a_ts(k, :), q_ts(k, :), ...
        ds(k), L, N_M, EPSILON);

    %% Finite difference solution
    ws_FD = load(sprintf("%s/FiniteDifference/composite/w_%d.mat", param_dir, k)).w_next;
    w_ts_FD = load(sprintf("%s/FiniteDifference/composite/w_t_%d.mat", param_dir, k)).w_t;
    ps_FD = load(sprintf("%s/FiniteDifference/composite/p_%d.mat", param_dir, k)).p;

    %% Plot displacement
    nexttile(1);
    plot(xs, -ws);
    hold on;
    plot(xs(1 : end - 1), -ws_FD);
    hold off;
    xlim([0, xMax]);

    %% Plot velocity
    nexttile(2);
    plot(xs, -w_ts);
    hold on;
    plot(xs(1 : end - 1), -w_ts_FD);
    hold off;
    xlim([0, xMax]);

    %% Plot pressure
    nexttile(3);
    plot(xs, ps);
    hold on;
    plot(xs(1 : end - 1), ps_FD);
    hold off;
    xlim([0, xMax]);
    
    %% Draw
    drawnow;
    pause(0.01);
end