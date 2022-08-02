clear;
close all;

% Adds analytical scripts to path
AnalyticalDir = "../Analytical_Scripts/Imposed/"
addpath(AnalyticalDir);


%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0,'defaultLegendFontSize', 18, 'DefaultLegendFontSizeMode','manual');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(groot, 'DefaultLegendInterpreter', 'latex');

% Main figure options
width = 900;
height = 450;

%% Load in color map
mapObj = load(sprintf("%s/Figure_Scripts/fine_red_blue_cmap.mat", AnalyticalDir));
cmap = mapObj.cmap;
blueCol = cmap(1, :);
redCol = cmap(end, :);
blackCol = [0, 0, 0];


%% Computational parameters
DELTA_T = 1e-4; % Timestep
IMPACT_TIME = 0.125; % Time of impact
T_MAX = 0.8; % Maximum time to plot
tsAnalytical = 0 : DELTA_T : T_MAX - IMPACT_TIME; % Analytical time array
MAX_TIMESTEP = T_MAX / DELTA_T; % Maximum timestep
NO_TIMESTEPS = length(tsAnalytical); % Number of timesteps
freq = 100; % Frequency of scatter
dimension = "2D";
type = "curvedDNS";
[epsilon, q, omega, p, L] = substrateparameters(type)


%%
fig = tiledlayout(1, 3);

% Need the flat functions struct
FlatFunctions = substratefunctions(type, dimension);

% Substrate position
nexttile;
ws = FlatFunctions.a(tsAnalytical);
plot(tsAnalytical, ws, 'color', blackCol, 'linewidth', 2);
xlabel("$t$");
ylabel("$w(t$)");
grid on;

% Substrate velocity
nexttile;
w_ts = FlatFunctions.a_t(tsAnalytical);
plot(tsAnalytical, w_ts, 'color', blackCol, 'linewidth', 2);
xlabel("$t$");
ylabel("$\dot{w}(t)$");
grid on;

% Substrate acceleration
nexttile;
w_tts = FlatFunctions.a_tt(tsAnalytical);
plot(tsAnalytical, w_tts, 'color', blackCol, 'linewidth', 2);
xlabel("$t$");
ylabel("$\ddot{w}(t)$");
grid on;

% Figure settings
set(gcf,'position', [100, 100, 800, 250]);
pause(0.1);


%%
figure(5);
hold on;
CurvedFunctions = substratefunctions(type, dimension);

xs = linspace(0, L, 1e3);

% Determine number of times
freq = 250;
tIdxs = 1 : freq : length(tsAnalytical);

colorFreq = floor(length(cmap) / length(tIdxs));

for k = 1 : length(tIdxs)
    t = tsAnalytical(tIdxs(k));
    if t == 0
        ws = zeros(size(xs));
    else
        ws = CurvedFunctions.w(xs, t);
    end
    plot(xs, ws, 'color', cmap(k * colorFreq, :), 'linewidth', 1.5);
end
grid on;
xlabel("$x$");
ylabel("$w(x, t)$");
ylim([-0.6, 0.4]);

set(gcf,'position', [100, 100, 800, 400]);