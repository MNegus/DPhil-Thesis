function PressureEvolutionPlot(dimension)
%%PressureEvolution2D.m
%   Script to make figures comparing pressure (composite) along the
%   substrate in time for the stationary, plate and quadratic substrate
%   cases.

    close all;

    addpath("../");
    addpath("../Pressures");

    painters = true;

    %% Figure options
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    set(0,'defaultAxesFontSize', 18);
    set(0, 'defaultAxesTickLabelInterpreter', 'latex');
    set(0,'defaultLegendFontSize', 18, 'DefaultLegendFontSizeMode','manual');
    set(groot, 'DefaultLegendInterpreter', 'latex');

    if painters
        set(0, 'defaultFigureRenderer', 'painters');
    else
        set(0, 'defaultFigureRenderer', 'opengl');
    end

    %% Load in color map
    mapObj = load("../../fine_red_blue_cmap.mat");
    cmap = mapObj.cmap;
    blueCol = cmap(1, :);
    redCol = cmap(end, :);

    %% Load parameters
    % Substrate parameters
    [epsilon, ~, ~, ~] = substrateparameters(); 

    tmax = 1;
    ts = linspace(0.05, tmax, 1000)';
    xMax = 1.25 * 2 * epsilon * sqrt(tmax);
    noPoints = 1e3;
    
    % Types of substrate
    if dimension == "2D"
        typeArr = ["stationary", "flat", "curved"];
        displayNames = ["Stationary substrate", ...
            "Moving substrate (rigid)", "Moving substrate (curved)"];
    else
        typeArr = ["stationary", "flat"];
        displayNames = ["Stationary substrate", ...
            "Moving substrate"];
    end

    %% Pressure in time plot
    freq = 100; % How frequent the lines are

    tiledlayout(length(typeArr), 1);

    tileNo = 1;
    for typeIdx = 1 : length(typeArr)
        type = typeArr(typeIdx);

        % Set line colors
        if type == "stationary"
                lineColor = 'black';

        elseif type == "flat"
            lineColor = redCol;

        else
            lineColor = blueCol;
        end

        % Load substrate functions
        SubstrateFunctions = imposedsubstratefunctions(type, dimension);

        % Move to next tile
        nexttile(tileNo);
        hold on;
        tileNo = tileNo + 1;

        for k = 1 : freq : length(ts)
            t = ts(k);

            %% Determine x values, clustered around the maximum pressure point
            % Find maximum pressure point
            [~, x_pMax] = pressuremax(t, SubstrateFunctions);

            fineWidth = 1e-3;
            leftWidth = x_pMax/ xMax;
            rightWidth = (xMax - x_pMax) / xMax;
            xsLower = linspace(0, x_pMax - fineWidth, floor(leftWidth * noPoints));
            xsFine = linspace(x_pMax - fineWidth, x_pMax + fineWidth, 100);

            xsUpper = linspace(x_pMax + fineWidth, xMax, floor(rightWidth * noPoints));

            xs = [xsLower, xsFine, xsUpper];

            % Load in composite pressure
            [ps, ~, ~] = substratepressure(xs, t, SubstrateFunctions);

            % Plot pressure
            h(1) = plot(xs, ps, 'Linewidth', 1.2, 'color', lineColor, ...
                'Displayname', displayNames(typeIdx));

        end

        %% Plot maximum pressure
        tLongs = linspace(0, 2 * tmax);
        [pMaxs, xMaxs] = pressuremax(tLongs, SubstrateFunctions);
        plot(xMaxs, pMaxs, 'color', lineColor, 'linestyle', '--');

        %% Tile figure settings
        xlim([0, xMax]);
        ylim([0, 500]);
        legend(h(1), 'location', 'northeast');
        if dimension == "2D"
            xlabel("$x$");
            ylabel("$p_{{comp}}(x, t)$");
            dirName = "Two-dimensional_Figures";
        else
            xlabel("$r$");
            ylabel("$p_{{comp}}(r, t)$");
            dirName = "Axisymmetric_Figures";
        end
        grid on;

    end

    %% Figure settings
    
    if dimension == "2D"
        dirName = "Two-dimensional_Figures";
        height = 600;
    else
        dirName = "Axisymmetric_Figures";
        height = 400;
    end
    
    set(gcf,'position', [100, 100, 900, height]);
    pause(0.1);
    
    % Export figure
    filename = sprintf("PressureEvolution_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);

end
