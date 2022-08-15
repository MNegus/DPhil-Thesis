%% AnalyticalParameterCompare.m
%
clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Forces/");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;

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
% tMax = 0.125;

% Define tMax to be where d = 0.25
tMax = (0.5)^2 / 3;

%% Load stationary solution
StatStruct = load("AnalyticalSolutions/StationarySol.mat").SolStruct;
dMaxStat = StatStruct.SubstrateFunctions.d(tMax);
d_tMaxStat = StatStruct.SubstrateFunctions.d_t(tMax);
HMaxStat = (1 + 4 / pi) * StatStruct.SubstrateFunctions.J(tMax);


%% Loops over alphas and compares turnovers
% ALPHAS = linspace(1, 200, 100);
% ALPHAS = 10.^linspace(-1, log10(200), 100);
ALPHAS = 10.^linspace(-2, log10(300), 500);
BETA = 0;
GAMMA = 0;

dMaxs = zeros(size(ALPHAS));
d_tMaxs = zeros(size(ALPHAS));
HMaxs = zeros(size(ALPHAS));

for ALPHAIdx = 1 : length(ALPHAS)
    ALPHA = ALPHAS(ALPHAIdx);

    % Load analytical solution
    fileName = append("AnalyticalSolutions/Fine_ALPHA_varying/ALPHA_", num2str(ALPHA), ".mat");
    SolStruct = load(fileName).SolStruct;

    % Save maximum ds and 
    dMaxs(ALPHAIdx) = SolStruct.SubstrateFunctions.d(tMax);
    d_tMaxs(ALPHAIdx) = SolStruct.SubstrateFunctions.d_t(tMax);
    HMaxs(ALPHAIdx) = (1 + 4 / pi) * SolStruct.SubstrateFunctions.J(tMax);
end

% Plot
tiledlayout(1, 3);

nexttile;
plot(ALPHAS, dMaxs, '-o');
set(gca, 'XScale', 'log');
grid on;
xlabel("$\alpha$");
ylabel("$d_{max}$");

nexttile;
plot(ALPHAS, d_tMaxs, '-o');
set(gca, 'XScale', 'log');
grid on;
xlabel("$\alpha$");
ylabel("$\dot{d}_{max}$");

nexttile;
plot(ALPHAS, HMaxs, '-o');
set(gca, 'XScale', 'log');
grid on;
xlabel("$\alpha$");
ylabel("$H_{max}$");

%% Loops over betas and compares turnovers
ALPHA = 2;
GAMMA = 100;
BETA_crit = 2 * sqrt(ALPHA * GAMMA);
% BETAS = BETA_crit * 10.^linspace(-2, 2, 100);
BETAS = BETA_crit * 10.^linspace(-3, 3, 500);

ts = linspace(0, tMax, 1e2);
dMaxs = zeros(size(BETAS));
d_tMaxs = zeros(size(BETAS));
HMaxs = zeros(size(BETAS));

for BETAIdx = 1 : length(BETAS)
    BETA = BETAS(BETAIdx);

    % Load analytical solution
    fileName = append("AnalyticalSolutions/Fine_BETA_varying/BETA_", num2str(BETA), ".mat");
    SolStruct = load(fileName).SolStruct;

    % Save maximum ds and 
    dMaxs(BETAIdx) = SolStruct.SubstrateFunctions.d(tMax);
    d_tMaxs(BETAIdx) = SolStruct.SubstrateFunctions.d_t(tMax);
    HMaxs(BETAIdx) = (1 + 4 / pi) * SolStruct.SubstrateFunctions.J(tMax);
end

% Plot
tiledlayout(1, 3);

nexttile;
plot(BETAS, dMaxs, '-o');
hold on;
xline(BETA_crit, 'LineStyle','--');
set(gca, 'XScale', 'log');
grid on;
xlabel("$\beta$");
ylabel("$d_{max}$");

nexttile;
plot(BETAS, d_tMaxs, '-o');
set(gca, 'XScale', 'log');
grid on;
xlabel("$\beta$");
ylabel("$\dot{d}_{max}$");

nexttile;
plot(BETAS, HMaxs, '-o');
set(gca, 'XScale', 'log');
grid on;
xlabel("$\beta$");
ylabel("$H_{max}$");

%% Loops over gammas and compares turnovers
ALPHA = 2;
BETA = 0;
% GAMMAS = 10.^linspace(-1, 6, 500);
GAMMAS = 10.^linspace(-1, 7, 500);

ts = linspace(0, tMax, 1e2);
dMaxs = zeros(size(GAMMAS));
d_tMaxs = zeros(size(GAMMAS));
HMaxs = zeros(size(GAMMAS));

for GAMMAIdx = 1 : length(GAMMAS)
    GAMMA = GAMMAS(GAMMAIdx);

    % Load analytical solution
    fileName = append("AnalyticalSolutions/Fine_GAMMA_varying/GAMMA_", num2str(GAMMA), ".mat");
    SolStruct = load(fileName).SolStruct;

    % Save maximum ds and 
    dMaxs(GAMMAIdx) = SolStruct.SubstrateFunctions.d(tMax);
    d_tMaxs(GAMMAIdx) = SolStruct.SubstrateFunctions.d_t(tMax);
    HMaxs(GAMMAIdx) = (1 + 4 / pi) * SolStruct.SubstrateFunctions.J(tMax);
end

% Plot
tiledlayout(1, 3);

nexttile;
plot(GAMMAS, dMaxs);
hold on;
yline(dMaxStat);

set(gca, 'XScale', 'log');
grid on;
xlabel("$\gamma$");
ylabel("$d_{max}$");

nexttile;
plot(GAMMAS, d_tMaxs);
hold on;
yline(d_tMaxStat);

set(gca, 'XScale', 'log');
grid on;
xlabel("$\gamma$");
ylabel("$\dot{d}_{max}$");

nexttile;
plot(GAMMAS, HMaxs);
hold on;
yline(HMaxStat);

set(gca, 'XScale', 'log');
grid on;
xlabel("$\gamma$");
ylabel("$H_{max}$");

%% Plot maximum gamma solution
GAMMA = GAMMAS(end);
fileName = append("AnalyticalSolutions/Fine_GAMMA_varying/GAMMA_", num2str(GAMMA), ".mat");
SolStruct = load(fileName).SolStruct;

figure(3);
plot(SolStruct.ts, SolStruct.ws, '-o');