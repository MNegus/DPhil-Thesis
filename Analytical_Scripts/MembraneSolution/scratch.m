%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD

addpath("FiniteDifference/");
addpath("FiniteDifference/PressuresFD/");
addpath("NormalModes/");

parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";


%% Parameters
EPSILON = 1;
% ALPHA = 2.004; BETA = 0 * EPSILON^2; GAMMA = 1069; 
delta = 1.0;
ALPHA = 1.1 * delta; BETA = 0.0; GAMMA = 668 * delta^3; 
L = 16;
T_MAX = 0.3;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 10924;

xs = linspace(0, L, N_MEMBRANE)';

ts_analytical = 0 : DELTA_T : T_MAX;

%% Finite differences
fd_data_dir = sprintf("%s/FiniteDifference", parent_dir);
SaveFDSolution(fd_data_dir, ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, "composite");


%% Validated NM
nm_data_dir = sprintf("%s/NormalModes", parent_dir);
SaveValidatedNMSolution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX, DELTA_T);

%% Load NM solutions
SolStruct = load(sprintf("%s/NormalModes/SolStruct.mat", parent_dir)).SolStruct;
as = SolStruct.as;
a_ts = SolStruct.a_ts;
q_ts = SolStruct.q_ts;
ds = SolStruct.ds;
N_M = SolStruct.N;


%% Plot in time
figure(1);

for k = 2 : 10 : length(ts_analytical)
    t = ts_analytical(k)
    %% Normal modes solution
    [ws, w_ts, ps] ...
        = MembraneSolutionNM(xs, as(k, :), ...
        a_ts(k, :), q_ts(k, :), ...
        ds(k), L, N_M, EPSILON);

    %% Finite difference solution
    ws_FD = load(sprintf("%s/FiniteDifference/composite/w_%d.mat", parent_dir, k)).w_next;

    %% Plot
    plot(xs, -ws);
    hold on;
    plot(xs(1 : end - 1), -ws_FD);
    hold off;
    drawnow;
    pause(0.01);
end


%% Plot turnovers
% Load FD ds
ds_FD = load(sprintf("%s/FiniteDifference/composite/ds.mat", parent_dir)).ds;

close(figure(2))
figure(2);
hold on;
plot(ts_analytical, 2 * sqrt(ts_analytical));
plot(ts_analytical, ds);
plot(ts_analytical, ds_FD);

%% Plot origin pressure

