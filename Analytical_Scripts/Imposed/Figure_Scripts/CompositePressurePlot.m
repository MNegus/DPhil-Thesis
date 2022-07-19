function CompositePressurePlot(dimension)
%%CompositePressure2D.m
%   Script to make a figure comparing the composite pressure to the outer
%   and inner pressures in 2D.

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

    %% Load parameters
    % Substrate parameters
    [epsilon, ~, ~, ~] = substrateparameters(); 

    xs = linspace(0, 0.4, 5e3);

    t = 2;

    %% Load in substrate function
    SubstrateFunctions = substratefunctions("stationary", dimension);

    %% Find pressures
    [ps_composite, ps_outer, ps_inner] ...
        = substratepressure(xs, t, SubstrateFunctions);

    % Restrict outer solution where zero
    ps_outer(ps_outer == 0) = nan;

    %% Set line colour
    lineColor = 'Black';

    %% Plot pressures
    figure(1);
    hold on;

    % Composite pressure
    h(3) = plot(xs, ps_composite, 'Linewidth', 2, 'Color', lineColor, ...
            'Displayname', 'Composite solution');

    % Outer pressure
    h(1) = plot(xs, ps_outer, 'Linewidth', 3, 'Color', redCol, ...
        'Linestyle', '--', 'Displayname', 'Outer solution');

    % Inner pressure
    h(2) = plot(xs, ps_inner, 'Linewidth', 3, 'Color', blueCol, ...
        'Linestyle', ':', 'Displayname', 'Inner solution');


    %% Plot turnover line
    xline(epsilon * SubstrateFunctions.d(t), 'Linewidth', 0.5, ...
        'Color', 0.5 * [1 1 1], 'Linestyle', '--');

    %% Figure options
    grid on;
    box on;

    % Axes labels
    if dimension == "2D"
        dirName = "Two-dimensional_Figures";
        xlabel("$x$");
        ylabel("$p(x, -\epsilon^2 w(x, t), t)$");
    else
        dirName = "Axisymmetric_Figures";
        xlabel("$r$");
        ylabel("$p(r, 0, t)$");
    end

    % Set axes limits
    xlim([0, 0.35]);
    ylim([0, 40]);

    % Legend
    legend(h(1:3), 'Location', 'Northwest');

    % Set figure position
    set(gcf,'position', [100, 100, 800, 400]);

    set(gcf, 'Renderer', 'Painters');
    pause(0.1);

    filename = sprintf("CompositePressure_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);

end

