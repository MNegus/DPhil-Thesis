%% save_solutions.m
% Given parameters, saves the solutions using normal modes and FD

addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/Pressures");
addpath("../Analytical_Scripts/MembraneSolution/FiniteDifference/");
addpath("../Analytical_Scripts/MembraneSolution/FiniteDifference/PressuresFD/");
addpath("../Analytical_Scripts/MembraneSolution/NormalModes/");
addpath("../Analytical_Scripts/ImposedSolution/")

parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";

%% Parameters
EPSILON = 1;
L = 16;
T_MAX = 0.35;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 21848;

DELTA_X = L / (N_MEMBRANE - 1); 
M = N_MEMBRANE - 1;
xs = (0 : DELTA_X : L - DELTA_X)';

ts_analytical = 0 : DELTA_T : T_MAX;

%% Stationary membrane case
stationary_dir = sprintf("%s/Stationary", parent_dir);
mkdir(stationary_dir);

fd_data_dir = sprintf("%s/%s", stationary_dir, "composite");
mkdir(fd_data_dir);

% Save substrate functions
SubstrateFunctions = imposedsubstratefunctions("stationary", "2D");
SubstrateFunctions.epsilon = EPSILON;
SubstrateFunctions.L = L;
SubstrateFunctions = substratedependents(SubstrateFunctions);


% Save initial solution
w_next = zeros(size(xs));
w_t = zeros(size(xs));

for k = 1 : length(ts_analytical)
    t = ts_analytical(k);
    t
    % Determine pressure
    if t == 0
        p = zeros(size(xs));
    else
        [p, ~, ~] = substratepressure(xs', t, SubstrateFunctions);
        p = p';
    end
    plot(xs, p)
    xlim([0, 1])
    drawnow

    % Save solutions
    save(sprintf("%s/w_%d.mat", fd_data_dir, k), 'w_next');
    save(sprintf("%s/w_t_%d.mat", fd_data_dir, k), 'w_t');
    save(sprintf("%s/p_%d.mat", fd_data_dir, k), 'p');
end


%% Rubber sheet case
delta = 1.0;
ALPHA = 1.1 * delta; BETA = 0; GAMMA = 668.0 * delta^3;
rubber_sheet_dir = sprintf("%s/RubberSheet", parent_dir);
mkdir(rubber_sheet_dir);

% Normal modes
nm_data_dir = sprintf("%s/NormalModes", rubber_sheet_dir);
mkdir(nm_data_dir);
SaveValidatedNMSolution(nm_data_dir, ALPHA, BETA, GAMMA, EPSILON, L, T_MAX, DELTA_T);

% Finite difference
fd_data_dir = sprintf("%s/FiniteDifference", rubber_sheet_dir);
mkdir(fd_data_dir);
mkdir(append(fd_data_dir, "/composite"));
SaveFDSolution(fd_data_dir, ALPHA, BETA, GAMMA, EPSILON, N_MEMBRANE, L, T_MAX, DELTA_T, "composite");


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