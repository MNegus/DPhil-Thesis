%% AnalyticalParameterCompare.m

clear;
close all;


addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/Pressures");
addpath("../Analytical_Scripts/MembraneSolution/FiniteDifference/");
addpath("../Analytical_Scripts/MembraneSolution/FiniteDifference/PressuresFD/");
addpath("../Analytical_Scripts/MembraneSolution/NormalModes/");
addpath("../Analytical_Scripts/ImposedSolution/")

parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";

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
EPSILON = 1;
L = 16;
T_MAX = 0.35;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 21848;

DELTA_X = L / (N_MEMBRANE - 1); 
M = N_MEMBRANE - 1;
xs = (0 : DELTA_X : L - DELTA_X)';

ts_analytical = 0 : DELTA_T : T_MAX

tTest = 1 / 16;
tTestIdx = tTest / DELTA_T;
% tTestIdx = find(ts_analytical == tTest)

%% Loop over GAMMAS
BETA = 0;
ALPHA = 1.1;
GAMMAS = 668.0 * 10.^linspace(-3, 2, 101);
gamma_vary_dir = sprintf("%s/GAMMA_varying", parent_dir);

GAMMA_TESTS = GAMMAS(10: 20 : end);

% Arrays for solutions
dsMax = zeros(size(GAMMAS));
d_tsMax = zeros(size(GAMMAS));
JsMax = zeros(size(GAMMAS));
wsMax = zeros(size(GAMMAS));


for GAMMAIdx = 1 : length(GAMMAS)

    GAMMA = GAMMAS(GAMMAIdx)

    % GAMMA directory
    param_dir = sprintf("%s/GAMMA_%g", gamma_vary_dir, GAMMA);

    % Normal modes
    nm_data_dir = sprintf("%s/NormalModes", param_dir);

    % Load solution struct
    SolStruct = load(sprintf("%s/SolStruct.mat", nm_data_dir)).SolStruct;
    as = SolStruct.as;
    a_ts = SolStruct.a_ts;
    q_ts = SolStruct.q_ts;
    N_M = SolStruct.N;
    ds = SolStruct.ds;
    d_ts = SolStruct.d_ts;


    % Save maximum variables
    dsMax(GAMMAIdx) = interp1(SolStruct.ts, SolStruct.ds, tTest);
    d_tsMax(GAMMAIdx) = interp1(SolStruct.ts, SolStruct.d_ts, tTest);
    JsMax(GAMMAIdx) = interp1(SolStruct.ts, SolStruct.Js, tTest);

    % Save maximum displacement
    [wsMax(GAMMAIdx), ~, ~] = MembraneSolutionNM(0, as(tTestIdx, :), ...
            a_ts(tTestIdx, :), q_ts(tTestIdx, :), ds(tTestIdx), ...
            L, N_M, EPSILON);

end
%
% Plot gammas
fig = tiledlayout(2, 2);

nexttile(1);
scatter(GAMMAS, dsMax);
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
set(gca, 'XScale', 'Log');
ylabel("$d_0(t_c)$");
xlabel("$\gamma$");
xlim("padded");
ylim("padded");


nexttile(2);
scatter(GAMMAS, d_tsMax);
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
set(gca, 'XScale', 'Log');
ylabel("$\dot{d}_0(t_c)$")
xlabel("$\gamma$");
xlim("padded");
ylim("padded");

nexttile(3);
scatter(GAMMAS, JsMax);
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
set(gca, 'XScale', 'Log');
ylabel("$J(t_c)$");
xlabel("$\gamma$");
xlim("padded");
ylim("padded");

nexttile(4);
scatter(GAMMAS, wsMax);
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
set(gca, 'XScale', 'Log');
ylabel("$w(0, t_c)$");
xlabel("$\gamma$");
xlim("padded");
ylim("padded");



%% DELTA varying
BETA = 0;
DELTAS = 2.^linspace(-3, 3, 101);

DELTA_TESTS = 2.^linspace(-3, 3, 5)
ALPHA_TESTS = 1.1 * DELTA_TESTS
GAMMA_TESTS = 668.0 * DELTA_TESTS.^3

ALPHAS = 1.1 * DELTAS;
GAMMAS = 668.0 * DELTAS.^3;
delta_vary_dir = sprintf("%s/DELTA_varying", parent_dir);

% Arrays for solutions
dsMax = zeros(size(DELTAS));
d_tsMax = zeros(size(DELTAS));
JsMax = zeros(size(DELTAS));
wsMax = zeros(size(DELTAS));

for DELTAIdx = 1 : length(DELTAS)
    DELTA = DELTAS(DELTAIdx);
    ALPHA = ALPHAS(DELTAIdx);
    GAMMA = GAMMAS(DELTAIdx);

    % DELTA directory
    param_dir = sprintf("%s/DELTA_%g", delta_vary_dir, DELTA)

    % Normal modes
    nm_data_dir = sprintf("%s/NormalModes", param_dir);

    % Load solution struct
    SolStruct = load(sprintf("%s/SolStruct.mat", nm_data_dir)).SolStruct;

    as = SolStruct.as;
    a_ts = SolStruct.a_ts;
    q_ts = SolStruct.q_ts;
    N_M = SolStruct.N;
    ds = SolStruct.ds;
    d_ts = SolStruct.d_ts;


    % Save maximum variables
    dsMax(DELTAIdx) = interp1(SolStruct.ts, SolStruct.ds, tTest);
    d_tsMax(DELTAIdx) = interp1(SolStruct.ts, SolStruct.d_ts, tTest);
    JsMax(DELTAIdx) = interp1(SolStruct.ts, SolStruct.Js, tTest);

    % Save maximum displacement
    [wsMax(DELTAIdx), ~, ~] = MembraneSolutionNM(0, as(tTestIdx, :), ...
            a_ts(tTestIdx, :), q_ts(tTestIdx, :), ds(tTestIdx), ...
            L, N_M, EPSILON);

end

%
% Plot deltas
fig = tiledlayout(2, 2);

nexttile(1);
scatter(DELTAS, dsMax);
for DELTA = DELTA_TESTS
    xline(DELTA);
end
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$d_0(t_c)$");
xlabel("$\delta$");


nexttile(2);
scatter(DELTAS, d_tsMax);
for DELTA = DELTA_TESTS
    xline(DELTA);
end
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$\dot{d}_0(t_c)$")

nexttile(3);
scatter(DELTAS, JsMax);
for DELTA = DELTA_TESTS
    xline(DELTA);
end
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$J(t_c)$")

nexttile(4);
scatter(DELTAS, wsMax);
for DELTA = DELTA_TESTS
    xline(DELTA);
end
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$w(0, t_c)$")

%% BETA varying
BETAS = 10.^linspace(-1, 5, 101);
BETA_TESTS = 10.^linspace(-1, 5, 5);
ALPHA = 1.1;
GAMMA = 668.0;
beta_vary_dir = sprintf("%s/BETA_varying", parent_dir);

% Arrays for solutions
dsMax = zeros(size(BETAS));
d_tsMax = zeros(size(BETAS));
JsMax = zeros(size(BETAS));
wsMax = zeros(size(BETAS));

for BETAIdx = 1 : length(BETAS)
    BETA = BETAS(BETAIdx);

    % BETA directory
    param_dir = sprintf("%s/BETA_%g", beta_vary_dir, BETA);

    % Normal modes
    nm_data_dir = sprintf("%s/NormalModes", param_dir);

    % Load solution struct
    SolStruct = load(sprintf("%s/SolStruct.mat", nm_data_dir)).SolStruct;

    as = SolStruct.as;
    a_ts = SolStruct.a_ts;
    q_ts = SolStruct.q_ts;
    N_M = SolStruct.N;
    ds = SolStruct.ds;
    d_ts = SolStruct.d_ts;


    % Save maximum variables
    dsMax(BETAIdx) = interp1(SolStruct.ts, SolStruct.ds, tTest);
    d_tsMax(BETAIdx) = interp1(SolStruct.ts, SolStruct.d_ts, tTest);
    JsMax(BETAIdx) = interp1(SolStruct.ts, SolStruct.Js, tTest);

    % Save maximum displacement
    [wsMax(BETAIdx), ~, ~] = MembraneSolutionNM(0, as(tTestIdx, :), ...
            a_ts(tTestIdx, :), q_ts(tTestIdx, :), ds(tTestIdx), ...
            L, N_M, EPSILON);
end

%
% Plot betas
fig = tiledlayout(2, 2);

nexttile(1);
scatter(BETAS, dsMax);
for BETA = BETA_TESTS
    xline(BETA);
end
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$d_0(t_c)$");
xlabel("$\beta$");


nexttile(2);
scatter(BETAS, d_tsMax);
for BETA = BETA_TESTS
    xline(BETA);
end
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$\dot{d}_0(t_c)$")
xlabel("$\beta$");

nexttile(3);
scatter(BETAS, JsMax);
for BETA = BETA_TESTS
    xline(BETA);
end
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$J(t_c)$")
xlabel("$\beta$");


nexttile(4);
scatter(BETAS, wsMax);
for BETA = BETA_TESTS
    xline(BETA);
end
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$w(0, t_c)$")
xlabel("$\beta$");


%% Plot large delta/gamma solutions
% param_dir = sprintf("%s/DELTA_%g", delta_vary_dir, max(DELTAS));
% SolStruct = load(sprintf("%s/NormalModes/SolStruct.mat", param_dir)).SolStruct;
% as = SolStruct.as;
% a_ts = SolStruct.a_ts;
% q_ts = SolStruct.q_ts;
% ds = SolStruct.ds;
% N_M = SolStruct.N;
% d_ts = SolStruct.d_ts;
% Js = SolStruct.Js;
% 
% ts_analytical = 0 : DELTA_T : T_MAX;
% 
% figure(45);
% for k = 1 : length(ts_analytical)
%     t = ts_analytical(k)
% 
%     % Normal modes solution
%     [ws, w_ts, ps] ...
%         = MembraneSolutionNM(xs, as(k, :), ...
%         a_ts(k, :), ...
%         q_ts(k, :), ...
%         ds(k), L, N_M, EPSILON);
% 
%     % PLot solution
%     plot(xs, -ws);
%     pause(0.01);
% 
% end