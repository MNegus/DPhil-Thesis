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
sz = 20;

% Tiled layout sizes
width = 6;
height = 5;

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

tTest = 1 / 16;
tTestIdx = tTest / DELTA_T;
% tTestIdx = find(ts_analytical == tTest)

%% Stationary plate values
dStat = 2 * sqrt(tTest);
d_tStat = 1 / sqrt(tTest);
JStat = pi * dStat / (8 * d_tStat^2);
HStat = (1 + 4 / pi) * JStat;
pMaxStat = d_tStat^2 / 2;
EStat = pi * tTest;

%% Loop over GAMMAS
% close all;
% BETA = 0;
% ALPHA = 1.1;
% GAMMAS = 668.0 * 10.^linspace(-3, 2, 101);
% gamma_vary_dir = sprintf("%s/GAMMA_varying", parent_dir);
% 
% GAMMA_TESTS = GAMMAS(10: 20 : end);
% 
% xTicks = [10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5];
% 
% % Arrays for solutions
% dsMax = zeros(size(GAMMAS));
% d_tsMax = zeros(size(GAMMAS));
% psMax = zeros(size(GAMMAS));
% HsMax = zeros(size(GAMMAS));
% wsMax = zeros(size(GAMMAS));
% E_outersMax = zeros(size(GAMMAS));
% E_jetsMax = zeros(size(GAMMAS));
% N_Ms = zeros(size(GAMMAS));
% 
% 
% for GAMMAIdx = 1 : length(GAMMAS)
% 
%     GAMMA = GAMMAS(GAMMAIdx);
% 
%     % GAMMA directory
%     param_dir = sprintf("%s/GAMMA_%g", gamma_vary_dir, GAMMA);
% 
%     % Normal modes
%     nm_data_dir = sprintf("%s/NormalModes", param_dir);
% 
%     % Load solution struct
%     SolStruct = load(sprintf("%s/SolStruct.mat", nm_data_dir)).SolStruct;
%     as = SolStruct.as;
%     a_ts = SolStruct.a_ts;
%     q_ts = SolStruct.q_ts;
%     N_M = SolStruct.N;
%     ds = SolStruct.ds;
%     d_ts = SolStruct.d_ts;
%     E_outers = SolStruct.E_outers;
%     E_jets = SolStruct.E_jets;
% 
%     % Save maximum variables
%     dsMax(GAMMAIdx) = interp1(SolStruct.ts, SolStruct.ds, tTest);
%     d_tsMax(GAMMAIdx) = interp1(SolStruct.ts, SolStruct.d_ts, tTest);
%     psMax(GAMMAIdx) = d_tsMax(GAMMAIdx)^2 / 2;
%     HsMax(GAMMAIdx) = (1 + 4 / pi) * interp1(SolStruct.ts, SolStruct.Js, tTest);
%     E_outersMax(GAMMAIdx) = interp1(SolStruct.ts, E_outers, tTest);
%     E_jetsMax(GAMMAIdx) = interp1(SolStruct.ts, E_jets, tTest);
%     N_Ms(GAMMAIdx) = N_M;
% 
%     % Save maximum displacement
%     [wsMax(GAMMAIdx), ~, ~] = MembraneSolutionNM(0, as(tTestIdx, :), ...
%             a_ts(tTestIdx, :), q_ts(tTestIdx, :), ds(tTestIdx), ...
%             L, N_M, EPSILON);
% 
% end
% 
% % Plot gammas
% fig = tiledlayout(2, 2);
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * width, 0.5 * height, width, height]);
% 
% nexttile(1);
% scatter(GAMMAS, dsMax, sz, blueCol);
% for GAMMA = GAMMA_TESTS
%     xline(GAMMA);
% end
% yline(dStat, 'LineStyle', '--');
% set(gca, 'XScale', 'Log');
% ylabel("$d_0(t_c)$");
% xlabel("$\gamma$");
% xlim([10^-0.5, 10^5.5]);
% xticks(xTicks);
% set(gca,'xminorgrid','off','yminorgrid','off');
% ylim("padded");
% grid on;
% box on;
% 
% nexttile(2);
% scatter(GAMMAS, HsMax, sz, blueCol);
% for GAMMA = GAMMA_TESTS
%     xline(GAMMA);
% end
% yline(HStat, 'LineStyle', '--');
% set(gca, 'XScale', 'Log');
% ylabel("$H(t_c)$");
% xlabel("$\gamma$");
% xlim([10^-0.5, 10^5.5]);
% xticks(xTicks);
% set(gca,'xminorgrid','off','yminorgrid','off');
% ylim("padded");
% grid on;
% box on;
% 
% nexttile(3);
% scatter(GAMMAS, psMax, sz, blueCol);
% for GAMMA = GAMMA_TESTS
%     xline(GAMMA);
% end
% yline(pMaxStat, 'LineStyle', '--');
% set(gca, 'XScale', 'Log');
% ylabel("max($p(x, t_c)$)")
% xlabel("$\gamma$");
% xlim([10^-0.5, 10^5.5]);
% xticks(xTicks);
% set(gca,'xminorgrid','off','yminorgrid','off');
% ylim("padded");
% grid on;
% box on;
% 
% nexttile(4);
% hold on;
% scatter(GAMMAS, E_outersMax, sz, blueCol);
% scatter(GAMMAS, E_jetsMax, sz, redCol);
% for GAMMA = GAMMA_TESTS
%     xline(GAMMA);
% end
% yline(EStat, 'LineStyle', '--');
% set(gca, 'XScale', 'Log');
% ylabel("$E_K(t_c)$");
% xlabel("$\gamma$");
% xlim([10^-0.5, 10^5.5]);
% xticks(xTicks);
% set(gca,'xminorgrid','off','yminorgrid','off');
% ylim([0.095, 0.2]);
% grid on;
% box on;
% legend(["$E_{K, outer}$", "$E_{K, jets}$"], 'Location', 'Northwest');
% 
% figname = "MembraneFigures/MembraneGAMMAVary";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
% 
% % Plot the normal modes numbers
% figure(2);
% NHeight = 2;
% NWidth = 4;
% set(gcf,'units', 'inches', ...
%     'position',[0.5 * NWidth, 0.5 * NHeight, NWidth, NHeight]);
% scatter(GAMMAS, N_Ms, sz, blueCol);
% for GAMMA = GAMMA_TESTS
%     xline(GAMMA);
% end
% set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
% grid on;
% box on;
% set(gca,'xminorgrid','off','yminorgrid','off');
% ylabel("$N_M$");
% xlabel("$\gamma$");
% xlim([10^-0.5, 10^5.5]);
% ylim([10^1, 10^3]);
% xticks(xTicks);
% 
% figname = "MembraneFigures/NormalModesGAMMA";
% exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);


%% Loop over NEW GAMMAS
close all;
BETA = 0;
DELTA = 0.125;
ALPHA = 1.1 * DELTA;
GAMMAS = 10.^linspace(0, 6, 101);
gamma_vary_dir = sprintf("%s/NEW_GAMMA_varying", parent_dir);

GAMMA_TESTS = GAMMAS(1 : 25 : end);

xTicks = [10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6];
xLim = [10^-0.5, 10^6.5];

% Arrays for solutions
dsMax = zeros(size(GAMMAS));
d_tsMax = zeros(size(GAMMAS));
psMax = zeros(size(GAMMAS));
HsMax = zeros(size(GAMMAS));
wsMax = zeros(size(GAMMAS));
E_outersMax = zeros(size(GAMMAS));
E_jetsMax = zeros(size(GAMMAS));
N_Ms = zeros(size(GAMMAS));


for GAMMAIdx = 1 : length(GAMMAS)

    GAMMA = GAMMAS(GAMMAIdx);

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
    E_outers = SolStruct.E_outers;
    E_jets = SolStruct.E_jets;

    % Save maximum variables
    dsMax(GAMMAIdx) = interp1(SolStruct.ts, SolStruct.ds, tTest);
    d_tsMax(GAMMAIdx) = interp1(SolStruct.ts, SolStruct.d_ts, tTest);
    psMax(GAMMAIdx) = d_tsMax(GAMMAIdx)^2 / 2;
    HsMax(GAMMAIdx) = (1 + 4 / pi) * interp1(SolStruct.ts, SolStruct.Js, tTest);
    E_outersMax(GAMMAIdx) = interp1(SolStruct.ts, E_outers, tTest);
    E_jetsMax(GAMMAIdx) = interp1(SolStruct.ts, E_jets, tTest);
    N_Ms(GAMMAIdx) = N_M;

    % Save maximum displacement
    [wsMax(GAMMAIdx), ~, ~] = MembraneSolutionNM(0, as(tTestIdx, :), ...
            a_ts(tTestIdx, :), q_ts(tTestIdx, :), ds(tTestIdx), ...
            L, N_M, EPSILON);

end

% Plot gammas
fig = tiledlayout(2, 2);
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

nexttile(1);
scatter(GAMMAS, dsMax, sz, blueCol);
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
yline(dStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
ylabel("$d_0(t_c)$");
xlabel("$\gamma$");
xlim(xLim);
xticks(xTicks);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;

nexttile(2);
scatter(GAMMAS, HsMax, sz, blueCol);
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
yline(HStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
ylabel("$H(t_c)$");
xlabel("$\gamma$");
xlim(xLim);
xticks(xTicks);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;

nexttile(3);
scatter(GAMMAS, psMax, sz, blueCol);
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
yline(pMaxStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
ylabel("max($p(x, t_c)$)")
xlabel("$\gamma$");
xlim(xLim);
xticks(xTicks);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;

nexttile(4);
hold on;
scatter(GAMMAS, E_outersMax, sz, blueCol);
scatter(GAMMAS, E_jetsMax, sz, redCol);
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
yline(EStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
ylabel("$E_K(t_c)$");
xlabel("$\gamma$");
xlim(xLim);
xticks(xTicks);
set(gca,'xminorgrid','off','yminorgrid','off');
% ylim([0.095, 0.2]);
grid on;
box on;
legend(["$E_{K, outer}$", "$E_{K, jets}$"], 'Location', 'Northwest');

figname = "MembraneFigures/MembraneGAMMAVary";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Plot the normal modes numbers (NEW gamma varying)
% Determine phase and group velocity
ks = pi * (2 * N_Ms - 1) / (2 * L);
vPhases = sqrt((BETA + GAMMAS .* ks.^2) / ALPHA);
vGroups = (BETA + 2 * GAMMAS .* ks.^2) ./ sqrt(ALPHA * (BETA + GAMMAS .* ks.^2));


figure(2);
NHeight = 2;
NWidth = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * NWidth, 0.5 * NHeight, NWidth, NHeight]);

% Normal modes plot
h(1) = scatter(GAMMAS, N_Ms, sz, blueCol);
hold on;
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
h(2) = plot(GAMMAS, N_Ms(1) ./ GAMMAS.^0.25, 'Linewidth', lineWidth, 'Color', ...
    redCol, 'LineStyle', '--');
hold off;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;
box on;
set(gca,'xminorgrid','off','yminorgrid','off');
ylabel("$N_M$");
xlabel("$\gamma$");
xlim(xLim);
ylim('padded');
xticks(xTicks);
legend(h(1:2), ["NM solution", "$N_M \sim \gamma^{-1/4}$"], ...
    'Location', 'Northeast');

figname = "MembraneFigures/NormalModesGAMMA";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

% Velocity plot
figure(3);
set(gcf,'units', 'inches', ...
    'position',[0.5 * NWidth, 0.5 * NHeight, NWidth, NHeight]);
h(1) = scatter(GAMMAS, vPhases, sz, blueCol);
hold on;
h(2) = scatter(GAMMAS, vGroups, sz, redCol);
for GAMMA = GAMMA_TESTS
    xline(GAMMA);
end
hold off;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;
box on;
set(gca,'xminorgrid','off','yminorgrid','off');
ylabel("$v(k_{N_M})$");
xlabel("$\gamma$");
xlim(xLim);
ylim('padded');
xticks(xTicks);
legend(h(1:2), ["$v_{phase}(k_{N_M})$", "$v_{group}(k_{N_M})$"], ...
    'Location', 'Northwest');

figname = "MembraneFigures/PhaseGroupVelocitiesGAMMA";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);


%% DELTA varying
close all;
BETA = 0;
DELTAS = 2.^linspace(-3, 3, 101);

DELTA_TESTS = 2.^linspace(-3, 3, 5);
ALPHA_TESTS = 1.1 * DELTA_TESTS;
GAMMA_TESTS = 668.0 * DELTA_TESTS.^3;

xTicks = [1/8, 1/4, 1/2, 1, 2, 4, 8];
xLim = [2.^(-3.5), 2^(3.5)];

ALPHAS = 1.1 * DELTAS;
GAMMAS = 668.0 * DELTAS.^3;
delta_vary_dir = sprintf("%s/DELTA_varying", parent_dir);

% Arrays for solutions
dsMax = zeros(size(DELTAS));
d_tsMax = zeros(size(DELTAS));
psMax = zeros(size(DELTAS));
HsMax = zeros(size(DELTAS));
wsMax = zeros(size(DELTAS));
E_outersMax = zeros(size(DELTAS));
E_jetsMax = zeros(size(DELTAS));
N_Ms = zeros(size(DELTAS));

for DELTAIdx = 1 : length(DELTAS)
    DELTA = DELTAS(DELTAIdx);
    ALPHA = ALPHAS(DELTAIdx);
    GAMMA = GAMMAS(DELTAIdx);

    % DELTA directory
    param_dir = sprintf("%s/DELTA_%g", delta_vary_dir, DELTA);

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
    E_outers = SolStruct.E_outers;
    E_jets = SolStruct.E_jets;


    % Save maximum variables
    dsMax(DELTAIdx) = interp1(SolStruct.ts, SolStruct.ds, tTest);
    d_tsMax(DELTAIdx) = interp1(SolStruct.ts, SolStruct.d_ts, tTest);
    psMax(DELTAIdx) = d_tsMax(DELTAIdx)^2 / 2;
    HsMax(DELTAIdx) = (1 + 4 / pi) * interp1(SolStruct.ts, SolStruct.Js, tTest);
    E_outersMax(DELTAIdx) = interp1(SolStruct.ts, E_outers, tTest);
    E_jetsMax(DELTAIdx) = interp1(SolStruct.ts, E_jets, tTest);
    N_Ms(DELTAIdx) = N_M;

    % Save maximum displacement
    [wsMax(DELTAIdx), ~, ~] = MembraneSolutionNM(0, as(tTestIdx, :), ...
            a_ts(tTestIdx, :), q_ts(tTestIdx, :), ds(tTestIdx), ...
            L, N_M, EPSILON);

end

%
% Plot deltas
fig = tiledlayout(2, 2);
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

nexttile(1);
scatter(DELTAS, dsMax, sz, blueCol);
for DELTA = DELTA_TESTS
    xline(DELTA);
end
yline(dStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
ylabel("$d_0(t_c)$");
xlabel("$\delta$");
xticks(xTicks);
set(gca,'xminorgrid','off','yminorgrid','off');
xlim(xLim);
ylim("padded");
grid on;
box on;

nexttile(2);
scatter(DELTAS, HsMax, sz, blueCol);
for DELTA = DELTA_TESTS
    xline(DELTA);
end
yline(HStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$H(t_c)$")
xlabel("$\delta$");
xticks(xTicks);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;


nexttile(3);
scatter(DELTAS, psMax, sz, blueCol);
for DELTA = DELTA_TESTS
    xline(DELTA);
end
yline(pMaxStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("max($p(x, t_c)$)")
xlabel("$\delta$");
xticks(xTicks);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;


nexttile(4);
% scatter(DELTAS, wsMax);
hold on;
scatter(DELTAS, E_outersMax, sz, blueCol);
scatter(DELTAS, E_jetsMax, sz, redCol);
for DELTA = DELTA_TESTS
    xline(DELTA);
end
yline(EStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$E_K(t_c)$");
xlabel("$\delta$");
xticks(xTicks);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;
legend(["$E_{K, outer}$", "$E_{K, jets}$"], 'Location', 'Southeast');


%
figname = "MembraneFigures/MembraneDELTAVary";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Plot the normal modes numbers (NEW delta varying)
% Determine phase and group velocity
ks = pi * (2 * N_Ms - 1) / (2 * L);
vPhases = sqrt((BETA + GAMMAS .* ks.^2) ./ ALPHAS);
vGroups = (BETA + 2 * GAMMAS .* ks.^2) ./ sqrt(ALPHAS .* (BETA + GAMMAS .* ks.^2));

figure(2);
NHeight = 2;
NWidth = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * NWidth, 0.5 * NHeight, NWidth, NHeight]);

% Normal modes plot
h(1) = scatter(DELTAS, N_Ms, sz, blueCol);
hold on;
for DELTA = DELTA_TESTS
    xline(DELTA);
end
h(2) = plot(DELTAS, N_Ms(1) * (DELTAS(1) ./ DELTAS).^0.5, 'Linewidth', lineWidth, 'Color', ...
    redCol, 'LineStyle', '--');
hold off;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;
box on;
set(gca,'xminorgrid','off','yminorgrid','off');
ylabel("$N_M$");
xlabel("$\delta$");
xlim(xLim);
ylim([10^0.75, 10^2.5]);
xticks(xTicks);
legend(h(1:2), ["NM solution", "$N_M \sim \delta^{-1/2}$"], ...
    'Location', 'Northeast');

figname = "MembraneFigures/NormalModesDELTA";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

% Velocity plot
figure(3);
set(gcf,'units', 'inches', ...
    'position',[0.5 * NWidth, 0.5 * NHeight, NWidth, NHeight]);
h(1) = scatter(DELTAS, vPhases, sz, blueCol);
hold on;
h(2) = scatter(DELTAS, vGroups, sz, redCol);
for DELTA = DELTA_TESTS
    xline(DELTA);
end
hold off;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;
box on;
set(gca,'xminorgrid','off','yminorgrid','off');
ylabel("$v(k_{N_M})$");
xlabel("$\delta$");
xlim(xLim);
ylim([50, 2000]);
xticks(xTicks);
legend(h(1:2), ["$v_{phase}(k_{N_M})$", "$v_{group}(k_{N_M})$"], ...
    'Location', 'Northwest');

figname = "MembraneFigures/PhaseGroupVelocitiesDELTA";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);


%% BETA varying
close all;
BETAS = 3 * 10.^linspace(-1, 4, 101);
BETA_TESTS = 3 * 10.^linspace(-1, 4, 5);
BETA_TESTS = BETA_TESTS(2 : end);
ALPHA = 0.1375;
GAMMA = 1.305;
beta_vary_dir = sprintf("%s/BETA_varying", parent_dir);

% Axes labels and limits
% xTicks = [10^-1, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5];
xTicks = [10^-1, 10^1, 10^3, 10^5];
xLim = [10^-1.5, 10^5.5];

% Arrays for solutions
dsMax = zeros(size(BETAS));
d_tsMax = zeros(size(BETAS));
psMax = zeros(size(BETAS));
HsMax = zeros(size(BETAS));
wsMax = zeros(size(BETAS));
E_outersMax = zeros(size(BETAS));
E_jetsMax = zeros(size(BETAS));
N_Ms = zeros(size(BETAS));

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
    E_outers = SolStruct.E_outers;
    E_jets = SolStruct.E_jets;

    % Save maximum variables
    dsMax(BETAIdx) = interp1(SolStruct.ts, SolStruct.ds, tTest);
    d_tsMax(BETAIdx) = interp1(SolStruct.ts, SolStruct.d_ts, tTest);
    psMax(BETAIdx) = d_tsMax(BETAIdx)^2 / 2;
    HsMax(BETAIdx) = (1 + 4 / pi) * interp1(SolStruct.ts, SolStruct.Js, tTest);
    E_outersMax(BETAIdx) = interp1(SolStruct.ts, E_outers, tTest);
    E_jetsMax(BETAIdx) = interp1(SolStruct.ts, E_jets, tTest);
    N_Ms(BETAIdx) = N_M;

    % Save maximum displacement
    [wsMax(BETAIdx), ~, ~] = MembraneSolutionNM(0, as(tTestIdx, :), ...
            a_ts(tTestIdx, :), q_ts(tTestIdx, :), ds(tTestIdx), ...
            L, N_M, EPSILON);
end

%
% Plot betas
fig = tiledlayout(2, 2);
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

nexttile(1);
scatter(BETAS, dsMax, sz, blueCol);
for BETA = BETA_TESTS
    xline(BETA);
end
yline(dStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
ylim("padded");
ylabel("$d_0(t_c)$");
xlabel("$\beta$");
xticks(xTicks);
xlim(xLim);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;

nexttile(2);
scatter(BETAS, HsMax, sz, blueCol);
for BETA = BETA_TESTS
    xline(BETA);
end
yline(HStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("$H(t_c)$")
xlabel("$\beta$");
xticks(xTicks);
xlim(xLim);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;

nexttile(3);
scatter(BETAS, psMax, sz, blueCol);
for BETA = BETA_TESTS
    xline(BETA);
end
yline(pMaxStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
xlim("padded");
ylim("padded");
ylabel("max($p(x, t_c)$)")
xlabel("$\beta$");
xticks(xTicks);
xlim(xLim);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;


nexttile(4);
% scatter(BETAS, wsMax);
hold on;
scatter(BETAS, E_outersMax, sz, blueCol);
scatter(BETAS, E_jetsMax, sz, redCol);
for BETA = BETA_TESTS
    xline(BETA);
end
yline(EStat, 'LineStyle', '--');
set(gca, 'XScale', 'Log');
xlim(xLim);
ylim("padded");
ylabel("$E_K(t_c)$")
xlabel("$\beta$");
xticks(xTicks);
set(gca,'xminorgrid','off','yminorgrid','off');
ylim("padded");
grid on;
box on;
legend(["$E_{K, outer}$", "$E_{K, jets}$"], 'Location', 'Northwest');

figname = "MembraneFigures/MembraneBETAVary";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Plot the normal modes numbers (NEW beta varying)
% Determine phase and group velocity
ks = pi * (2 * N_Ms - 1) / (2 * L);
vPhases = sqrt((BETAS + GAMMA * ks.^2) / ALPHA);
vGroups = (BETAS + 2 * GAMMA * ks.^2) ./ sqrt(ALPHA * (BETAS + GAMMA .* ks.^2));

figure(2);
NHeight = 2;
NWidth = 3;
set(gcf,'units', 'inches', ...
    'position',[0.5 * NWidth, 0.5 * NHeight, NWidth, NHeight]);

% Normal modes plot
h(1) = scatter(BETAS, N_Ms, sz, blueCol);
hold on;
for BETA = BETA_TESTS
    xline(BETA);
end
h(2) = plot(BETAS, N_Ms(end - 3) * (BETAS(end - 3) ./ BETAS).^0.5, ...
    'Linewidth', lineWidth, 'Color', redCol, 'LineStyle', '--');
hold off;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;
box on;
set(gca,'xminorgrid','off','yminorgrid','off');
ylabel("$N_M$");
xlabel("$\beta$");
xlim(xLim);
ylim([30, 150]);
xticks(xTicks);
legend(h(1:2), ["NM solution", "$N_M \sim \beta^{-1/2}$"], ...
    'Location', 'Southwest');

figname = "MembraneFigures/NormalModesBETA";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

% Velocity plot
figure(3);
set(gcf,'units', 'inches', ...
    'position',[0.5 * NWidth, 0.5 * NHeight, NWidth, NHeight]);
h(1) = scatter(BETAS, vPhases, sz, blueCol);
hold on;
h(2) = scatter(BETAS, vGroups, sz, redCol);
for BETA = BETA_TESTS
    xline(BETA);
end
hold off;
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;
box on;
set(gca,'xminorgrid','off','yminorgrid','off');
ylabel("$v(k_{N_M})$");
xlabel("$\beta$");
xlim(xLim);
ylim('padded');
xticks(xTicks);
legend(h(1:2), ["$v_{phase}(k_{N_M})$", "$v_{group}(k_{N_M})$"], ...
    'Location', 'Northwest');

figname = "MembraneFigures/PhaseGroupVelocitiesBETA";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
