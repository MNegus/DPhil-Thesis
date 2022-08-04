close all;
clear;

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
mapObj = load("../fine_red_blue_cmap.mat");
cmap = mapObj.cmap;
blueCol = cmap(1, :);
redCol = cmap(end, :);
blackCol = [0, 0, 0];

%% Parameters
epsilon = 1;
ALPHAS = [1, 2];
BETAS = 0;
% GAMMAS = [1, 10, 100, 1000];
GAMMAS = [20];
tMax = 0.625;
ts = linspace(0, tMax, 1e3);
forceType = "composite";

%% Save stationary functions
StationarySubstrateFunctions ...
    = platesubstratefunctions(ts, zeros(size(ts)), zeros(size(ts)), ...
    zeros(size(ts)), epsilon);

%% Plot stationary solutions
tiledlayout(2, 3);

nexttile(4);
hold on;
plot(ts, StationarySubstrateFunctions.d(ts), 'linewidth', 2, ...
    'color', 'black', 'linestyle', '--');
xlabel("$t$");
ylabel("$d(t)$");
set(gca,'ColorOrderIndex',1);

nexttile(5);
hold on;
plot(ts, StationarySubstrateFunctions.d_t(ts), 'linewidth', 2, ...
    'color', 'black', 'linestyle', '--');
xlabel("$t$");
ylabel("$\dot{d}(t)$");
set(gca,'ColorOrderIndex',1);

nexttile(6);
hold on;
[force, ~, ~] = substrateforce(ts, StationarySubstrateFunctions);
plot(ts, force, 'linewidth', 2, ...
    'color', 'black', 'linestyle', '--');
xlabel("$t$");
ylabel("$F(t)$");
set(gca,'ColorOrderIndex',1);
drawnow;

%% Loop over params
for ALPHA = ALPHAS
    for BETA = BETAS
        for GAMMA = GAMMAS
            %% Solve ODE
            [ts, ws, w_ts, w_tts] = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, forceType);

            %% Plot solutions
            nexttile(1);
            hold on;
            plot(ts, ws, 'linewidth', 2);
            xlabel("$t$");
            ylabel("$w(t)$");
            
            nexttile(2);
            hold on;
            plot(ts, w_ts, 'linewidth', 2);
            xlabel("$t$");
            ylabel("$\dot{w}(t)$");
            
            nexttile(3);
            hold on;
            plot(ts, w_tts, 'linewidth', 2);
            xlabel("$t$");
            ylabel("$\dot{w}(t)$");

            %% Make substrate functions with solution
            SubstrateFunctions = platesubstratefunctions(ts, ws, w_ts, w_tts, epsilon);
            
            %% Plot turnover point
            nexttile(4);
            hold on;
            plot(ts, SubstrateFunctions.d(ts), 'linewidth', 2);
            
            %% Plot jet thickness
            nexttile(5);
            hold on;
            plot(ts, SubstrateFunctions.d_t(ts), 'linewidth', 2);
            ylim([0, 10]);
            
            %% Plot force 
            [force, ~, ~] = substrateforce(ts, SubstrateFunctions);
            nexttile(6);
            hold on;
            plot(ts, force, 'linewidth', 2);
            
            % Fig options
            set(gcf, 'position', [400, 400, 1000, 600]);
            drawnow;
        end
        
        
    end
end


%% TODO
% * Make outer and inner free surface plots to see the difference between
% the cases