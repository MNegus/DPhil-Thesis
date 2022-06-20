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


%% Parameters
xs = linspace(0, 0.275, 5e3);
epsilon = 0.1;

ts = [0.0625, 1.0];

%% Loop over the times and make the figures
figure(1);
hold on;

for tIdx = 1 : length(ts)
    t = ts(tIdx);
    
    %% Save structure arrays
    PlatePositions = platepositions(0, 0, 0);
    TimeDependents = timedependents(t, PlatePositions);
    
    %% Find pressures
    [ps_composite, ps_outer, ps_inner] ...
        = substratepressure(xs, PlatePositions, TimeDependents, epsilon);
    
    % Restrict outer solution where zero
    ps_outer(ps_outer == 0) = nan;
    
    %% Set line colour
    if tIdx == 2
        lineColor = 0.5 * [1 1 1];
    else
        lineColor = 'Black';
    end
    
    %% Plot pressures
    h(3 * (tIdx - 1) + 1) = plot(xs, ps_outer, 'Linewidth', 3, 'Color', lineColor, ...
        'Linestyle', '--', 'Displayname', 'Outer');
    h(3 * (tIdx - 1) + 2) = plot(xs, ps_inner, 'Linewidth', 3, 'Color', lineColor, ...
        'Linestyle', ':', 'Displayname', 'Inner');
    h(3 * (tIdx - 1) + 3) = plot(xs, ps_composite, 'Linewidth', 1.85, 'Color', lineColor, ...
        'Displayname', 'Composite');
    
    %% Plot turnover line
%     xline(epsilon * TimeDependents.ds, 'Linewidth', 0.5, 'Color', 0.5 * [1 1 1], ...
%         'Linestyle', '--');
    
end

%% Add annotations
annotation('textbox',...
    [0.14 0.5815 0.094 0.0675],...
    'String',{'$t = 0.0625$'}, 'interpreter', 'latex', 'Fontsize', 18, ...
    'Edgecolor', 'None');

annotation('textbox',...
    [0.65 0.5815 0.094 0.0675],...
    'String',{'$t = 1$'}, 'interpreter', 'latex', 'Fontsize', 18, ...
    'Edgecolor', 'None');

%% Figure options
grid on;
box on;

% Axes labels
xlabel("$x$");
ylabel("$p(x, 0, t)$");

% Set axes limits
xlim([0, 0.26]);
ylim([0, 200]);
xticks([0, 0.05, 0.1, 0.15, 0.2, 0.25]);

% Legend
legend(h(1:3), 'Location', 'North');


% Set figure position
set(gcf,'position', [100, 100, 960, 400]);

set(gcf, 'Renderer', 'Painters');

exportgraphics(gca,'png/CompositePressure2D.png', 'Resolution', 300);
exportgraphics(gca,'eps/CompositePressure2D.eps', 'Resolution', 300);
