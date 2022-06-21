%%CompositePressure2D.m
%   Script to make a figure comparing the composite pressure to the outer
%   and inner pressures in 2D.

clear;
close all;

addpath("../");
addpath("../Pressures");

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Load in color map
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;
blueCol = cmap(1, :);
redCol = cmap(end, :);

%% Parameters
xs = linspace(0, 0.4, 5e3);
epsilon = 0.1;

ts = [2];

%% Substrate definition
% Oscillatory time dependence
L = 1;
q = 0.1;
omega = 2.5;

a = @(t) q * L^2 * (t^2 + 1 - cos(omega * t));
a_t = @(t) q * L^2 * (2 * t + omega * sin(omega * t));
a_tt = @(t) q * L^2 * (2 + omega^2 * cos(omega * t));

b = @(t) - a(t) / L^2;
b_t = @(t) - a_t(t) / L^2;
b_tt = @(t) - a_tt(t) / L^2;


%% Loop over the times and make the figures
figure(1);
hold on;

for tIdx = 1 : length(ts)
    t = ts(tIdx);
    
    %% Save structure arrays
    % Stationary substrate case
    StationarySubstrateCoefficients = substratecoefficients(0, 0, 0, 0, 0, 0, epsilon);
    StationaryTimeDependents = timedependents(t, StationarySubstrateCoefficients);
    
    % Quadratic substrate case
    SubstrateCoefficients = substratecoefficients(a(t), b(t), a_t(t), ...
        b_t(t), a_tt(t), b_tt(t), epsilon);
    TimeDependents = timedependents(t, SubstrateCoefficients);
    
    %% Find stationary pressures
    [Stationary_ps_composite, Stationary_ps_outer, Stationary_ps_inner] ...
        = substratepressure(xs, StationarySubstrateCoefficients, StationaryTimeDependents, epsilon);
    
    % Restrict outer solution where zero
    Stationary_ps_outer(Stationary_ps_outer == 0) = nan;
    
    %% Set line colour

    StatlineColor = 'Black';
    lineColor = 0.5 * [1 1 1];
    
    %% Plot pressures
    
    % Composite pressures
    h(3) = plot(xs, Stationary_ps_composite, 'Linewidth', 2, 'Color', StatlineColor, ...
            'Displayname', 'Composite solution');
        
    % Outer pressures
    h(1) = plot(xs, Stationary_ps_outer, 'Linewidth', 3, 'Color', redCol, ...
        'Linestyle', '--', 'Displayname', 'Outer solution');
    
    % Inner pressures
    h(2) = plot(xs, Stationary_ps_inner, 'Linewidth', 3, 'Color', blueCol, ...
        'Linestyle', ':', 'Displayname', 'Inner solution');
    
    
%     
    %% Plot pressures
%     h(3 * (tIdx - 1) + 1) = plot(xs, Stationary_ps_outer, 'Linewidth', 3, 'Color', lineColor, ...
%         'Linestyle', '--', 'Displayname', 'Outer');
%     h(3 * (tIdx - 1) + 2) = plot(xs, Stationary_ps_inner, 'Linewidth', 3, 'Color', lineColor, ...
%         'Linestyle', ':', 'Displayname', 'Inner');
%     h(3 * (tIdx - 1) + 3) = plot(xs, Stationary_ps_composite, 'Linewidth', 1.85, 'Color', lineColor, ...
%         'Displayname', 'Composite');
    
    %% Plot turnover line
    xline(epsilon * StationaryTimeDependents.ds, 'Linewidth', 0.5, 'Color', 0.5 * [1 1 1], ...
        'Linestyle', '--');
    
end

%% Add annotations
% annotation('textbox',...
%     [0.14 0.5815 0.094 0.0675],...
%     'String',{'$t = 0.0625$'}, 'interpreter', 'latex', 'Fontsize', 18, ...
%     'Edgecolor', 'None');
% 
% annotation('textbox',...
%     [0.65 0.5815 0.094 0.0675],...
%     'String',{'$t = 1$'}, 'interpreter', 'latex', 'Fontsize', 18, ...
%     'Edgecolor', 'None');

%% Figure options
grid on;
box on;

% Axes labels
xlabel("$x$");
ylabel("$p(x, 0, t)$");

% Set axes limits
xlim([0, 0.35]);
ylim([0, 40]);
% xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25]);

% Legend
legend(h(1:3), 'Location', 'Northwest');


% Set figure position
set(gcf,'position', [100, 100, 600, 400]);

set(gcf, 'Renderer', 'Painters');
pause(0.1);
% 
exportgraphics(gca,'png/CompositePressure2D.png', 'Resolution', 300);
exportgraphics(gca,'eps/CompositePressure2D.eps', 'Resolution', 300);
