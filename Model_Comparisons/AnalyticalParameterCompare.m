%% AnalyticalParameterCompare.m
%
clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Forces/");
addpath("../Analytical_Scripts/Energies/")

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

% Scatter size
sz = 5;

% Tiled layout sizes
width = 6;
height = 5;


%% Parameters
% tMax = 0.125;

% Define tMax to be where d = 0.5
tMax = (0.5)^2 / 3;

%% Load stationary solution
StatStruct = load("AnalyticalSolutions/StationarySol.mat").SolStruct;
dMaxStat = StatStruct.SubstrateFunctions.d(tMax);
d_tMaxStat = StatStruct.SubstrateFunctions.d_t(tMax);
HMaxStat = (1 + 4 / pi) * StatStruct.SubstrateFunctions.J(tMax);
EOuterStat = interp1(StatStruct.ts, StatStruct.EOuters, tMax);
EJetStat = interp1(StatStruct.ts, StatStruct.EJets, tMax);

%% Loops over alphas and compares turnovers
% ALPHAS = linspace(1, 200, 100);
% ALPHAS = 10.^linspace(-1, log10(200), 100);
ALPHAS = 10.^linspace(-2, log10(300), 500);
BETA = 0;
GAMMA = 0;

ALPHA_TESTS = 2 * 10.^linspace(0, log10(100 / 2), 5);

dMaxs = zeros(size(ALPHAS));
d_tMaxs = zeros(size(ALPHAS));
HMaxs = zeros(size(ALPHAS));
EOuterMaxs = zeros(size(ALPHAS));
EJetMaxs = zeros(size(ALPHAS));
wMaxs = zeros(size(ALPHAS));

for ALPHAIdx = 1 : length(ALPHAS)
    ALPHA = ALPHAS(ALPHAIdx)

    % Load analytical solution
    fileName = append("AnalyticalSolutions/Fine_ALPHA_varying/ALPHA_", num2str(ALPHA), ".mat");
    SolStruct = load(fileName).SolStruct;

    % Save maximum ds and Hs
    dMaxs(ALPHAIdx) = SolStruct.SubstrateFunctions.d(tMax);
    d_tMaxs(ALPHAIdx) = SolStruct.SubstrateFunctions.d_t(tMax);
    HMaxs(ALPHAIdx) = (1 + 4 / pi) * SolStruct.SubstrateFunctions.J(tMax);

    % Save maximum energies
    EOuterMaxs(ALPHAIdx) = interp1(SolStruct.ts, SolStruct.EOuters, tMax);
    EJetMaxs(ALPHAIdx) = interp1(SolStruct.ts, SolStruct.EJets, tMax);
    
    wMaxs(ALPHAIdx) = SolStruct.SubstrateFunctions.w_t(tMax);
end

% Plot
tiledlayout(2, 2);
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

xTicks = 10.^(-2 : 2 : 2);

nexttile;
scatter(ALPHAS, dMaxs, sz, blueCol);
hold on;
yline(dMaxStat, 'LineStyle', '--');
for ALPHA = ALPHA_TESTS
    xline(ALPHA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-3, 10^3]);
xticks(xTicks);
ylim('padded');
xlabel("$\alpha$");
ylabel("$d(t_c)$");


nexttile;
scatter(ALPHAS, d_tMaxs, sz, blueCol);
hold on;
yline(d_tMaxStat, 'LineStyle', '--');
for ALPHA = ALPHA_TESTS
    xline(ALPHA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-3, 10^3]);
xticks(xTicks);
ylim('padded');
xlabel("$\alpha$");
ylabel("$\dot{d}(t_c)$");


nexttile;
scatter(ALPHAS, HMaxs, sz, blueCol);
hold on;
yline(HMaxStat, 'LineStyle', '--');
for ALPHA = ALPHA_TESTS
    xline(ALPHA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-3, 10^3]);
xticks(xTicks);
ylim('padded');
xlabel("$\alpha$");
ylabel("$H(t_c)$");


nexttile;
hold on;
h(1) = scatter(ALPHAS, EOuterMaxs, sz, blueCol);
h(2) = scatter(ALPHAS, EJetMaxs, sz, redCol);
yline(EOuterStat, 'LineStyle', '--');
for ALPHA = ALPHA_TESTS
    xline(ALPHA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-3, 10^3]);
xticks(xTicks);
ylim('padded');
xlabel("$\alpha$");
ylabel("$E(t_c)$");
legend(["$E_{K, outer}$", "$E_{K, splash}$"], 'Location', 'Southeast');

figname = "PlateFigures/ALPHAVaryAnalytical";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

% figure(7);
% plot(ALPHAS, wMaxs)
% set(gca, 'XScale', 'log')

%% Loops over betas and compares turnovers
ALPHA = 2;
GAMMA = 100;
BETA_crit = 2 * sqrt(ALPHA * GAMMA);
BETAS = BETA_crit * 10.^linspace(-3, 3, 500);

BETA_TESTS = BETA_crit * 10.^([-2, -1, 0, 1, 2]);
BETA_TESTS = [0.0, 7.07, 14.14, 28.28, 56.57, 113.1];

dMaxs = zeros(size(BETAS));
d_tMaxs = zeros(size(BETAS));
HMaxs = zeros(size(BETAS));
EOuterMaxs = zeros(size(BETAS));
EJetMaxs = zeros(size(BETAS));

xLimMax = 10^5

for BETAIdx = 1 : length(BETAS)
    BETA = BETAS(BETAIdx)

    % Load analytical solution
    fileName = append("AnalyticalSolutions/Fine_BETA_varying/BETA_", num2str(BETA), ".mat");
    SolStruct = load(fileName).SolStruct;

    % Save maximum ds and Hs
    dMaxs(BETAIdx) = SolStruct.SubstrateFunctions.d(tMax);
    d_tMaxs(BETAIdx) = SolStruct.SubstrateFunctions.d_t(tMax);
    HMaxs(BETAIdx) = (1 + 4 / pi) * SolStruct.SubstrateFunctions.J(tMax);

    % Save maximum energies
    EOuterMaxs(BETAIdx) = interp1(SolStruct.ts, SolStruct.EOuters, tMax);
    EJetMaxs(BETAIdx) = interp1(SolStruct.ts, SolStruct.EJets, tMax);
end

% Plot
tiledlayout(2, 2);
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, 1.1 * height]);

xTicks = 10.^(-2 : 2 : 5);

nexttile;
scatter(BETAS, dMaxs, sz, blueCol);
hold on;
yline(dMaxStat, 'LineStyle', '--');
for BETA = BETA_TESTS
    xline(BETA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-2, xLimMax]);
xticks(xTicks);
ylim([0.492, 0.501]);
xlabel("$\beta$");
ylabel("$d(t_c)$");


nexttile;
scatter(BETAS, d_tMaxs, sz, blueCol);
hold on;
yline(d_tMaxStat, 'LineStyle', '--');
for BETA = BETA_TESTS
    xline(BETA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-2, xLimMax]);
xticks(xTicks);
ylim([2.82, 3.05]);
xlabel("$\beta$");
ylabel("$\dot{d}(t_c)$");
% xticks(xTicks);


nexttile;
scatter(BETAS, HMaxs, sz, blueCol);
hold on;
yline(HMaxStat, 'LineStyle', '--');
for BETA = BETA_TESTS
    xline(BETA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-2, xLimMax]);
xticks(xTicks);
ylim([0.0192, 0.0202]);
xlabel("$\beta$");
ylabel("$H(t_c)$");
xticks(xTicks);


nexttile;
hold on;
h(1) = scatter(BETAS, EOuterMaxs, sz, blueCol);
h(2) = scatter(BETAS, EJetMaxs, sz, redCol);
yline(EOuterStat, 'LineStyle', '--');
for BETA = BETA_TESTS
    xline(BETA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-2, xLimMax]);
xticks(xTicks);
ylim([0.069, 0.085]);
xlabel("$\beta$");
ylabel("$E(t_c)$");
xticks(xTicks);

legend(["$E_{K, outer}$", "$E_{K, splash}$"], 'Location', 'best');

pause(0.5);
figname = "PlateFigures/BETAVaryAnalytical";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Loops over gammas and compares turnovers
ALPHA = 2;
BETA = 0;
GAMMAS = 10.^linspace(-1, 7, 500);

GAMMA_guaranteed = [9000, 21000];
a = log10(9000);
b = log10(21000);
inc = b - a;
GAMMA_TESTS = 10.^(a - 4 * inc : inc : a + 2 * inc)
GAMMA_TESTS(end)
%%
ts = linspace(0, tMax, 1e2);
dMaxs = zeros(size(GAMMAS));
d_tMaxs = zeros(size(GAMMAS));
HMaxs = zeros(size(GAMMAS));
EOuterMaxs = zeros(size(GAMMAS));
EJetMaxs = zeros(size(GAMMAS));

for GAMMAIdx = 1 : length(GAMMAS)
    GAMMA = GAMMAS(GAMMAIdx)

    % Load analytical solution
    fileName = append("AnalyticalSolutions/Fine_GAMMA_varying/GAMMA_", num2str(GAMMA), ".mat");
    SolStruct = load(fileName).SolStruct;

    % Save maximum ds and Hs
    dMaxs(GAMMAIdx) = SolStruct.SubstrateFunctions.d(tMax);
    d_tMaxs(GAMMAIdx) = SolStruct.SubstrateFunctions.d_t(tMax);
    HMaxs(GAMMAIdx) = (1 + 4 / pi) * SolStruct.SubstrateFunctions.J(tMax);

    % Save maximum energies
    EOuterMaxs(GAMMAIdx) = interp1(SolStruct.ts, SolStruct.EOuters, tMax);
    EJetMaxs(GAMMAIdx) = interp1(SolStruct.ts, SolStruct.EJets, tMax);
end

% Plot
tileFig = tiledlayout(2, 2);
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

xTicks = 10.^(-2 : 2 : 8);

nexttile;
scatter(GAMMAS, dMaxs, sz, blueCol);
hold on;
yline(dMaxStat, 'LineStyle', '--');
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-2, 10^8]);
xticks(xTicks);
xtickangle(0)
ylim([0.492, 0.501]);
xlabel("$\gamma$");
ylabel("$d(t_c)$");


nexttile;
scatter(GAMMAS, d_tMaxs, sz, blueCol);
hold on;
yline(d_tMaxStat, 'LineStyle', '--');
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-2, 10^8]);
xticks(xTicks);
xtickangle(0)
ylim([2.82, 3.05]);
xlabel("$\gamma$");
ylabel("$\dot{d}(t_c)$");


nexttile;
scatter(GAMMAS, HMaxs, sz, blueCol);
hold on;
yline(HMaxStat, 'LineStyle', '--');
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end

box on;
xlim([10^-2, 10^8]);
xticks(xTicks);
xtickangle(0)
ylim([0.0192, 0.0202]);


set(gca, 'XScale', 'log');
grid on;
xlabel("$\gamma$");
ylabel("$H(t_c)$");


nexttile;
hold on;
h(1) = scatter(GAMMAS, EOuterMaxs, sz, blueCol);
h(2) = scatter(GAMMAS, EJetMaxs, sz, redCol);
yline(EOuterStat, 'LineStyle', '--');
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end

set(gca, 'XScale', 'log');
grid on;
box on;
xlim([10^-2, 10^8]);
xticks(xTicks);
xtickangle(0)
ylim([0.069, 0.085]);
xlabel("$\gamma$");
ylabel("$E(t_c)$");
legend(["$E_{K, outer}$", "$E_{K, splash}$"], 'Location', 'Northwest');


figname = "PlateFigures/GAMMAVaryAnalytical";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Plot maximum gamma solution
% figure(3);
% hold on;
% plot(SolStruct.ts, sqrt(3) ./ (2 * sqrt(SolStruct.ts)))
% 
% 
% GAMMA = GAMMAS(242)
% fileName = append("AnalyticalSolutions/Fine_GAMMA_varying/GAMMA_", num2str(GAMMA), ".mat");
% SolStruct = load(fileName).SolStruct;
% plot(SolStruct.ts, SolStruct.d_ts);
% 
% 
% GAMMA = GAMMAS(331)
% fileName = append("AnalyticalSolutions/Fine_GAMMA_varying/GAMMA_", num2str(GAMMA), ".mat");
% SolStruct = load(fileName).SolStruct;
% plot(SolStruct.ts, SolStruct.d_ts);

%% Plot energies
% figure(4);
% hold on;
% GAMMA_energ = [GAMMAS(276), GAMMAS(310), GAMMAS(395)];
% 
% % for GAMMAIdx = 1 : length(GAMMA_energ)
% for GAMMAIdx = 1
%     GAMMA = GAMMA_energ(GAMMAIdx)
% 
%     % Load analytical solution
%     fileName = append("AnalyticalSolutions/Fine_GAMMA_varying/GAMMA_", num2str(GAMMA), ".mat");
%     SolStruct = load(fileName).SolStruct;
%     ts = SolStruct.ts;
% 
%     % Plot maximum energy
%     plot(SolStruct.ts, SolStruct.EOuters)
%     plot(SolStruct.ts, SolStruct.EJets, 'LineStyle', '--');
% end
% 
% plot(ts, 2 * sqrt(3) * ts.^(3/2), 'color', 'black', 'LineStyle', '--', ...
%     'LineWidth', 2);

