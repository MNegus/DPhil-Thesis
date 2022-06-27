function JetEvolutionPlot(dimension)
%JETEVOLUTIONPLOT
% Script to plot the time evolution of the jets/splash sheet

    close all;

    addpath("../");

    %% Figure options
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    set(0,'defaultAxesFontSize', 18);
    set(0, 'defaultAxesTickLabelInterpreter', 'latex');
    set(0, 'defaultFigureRenderer', 'painters');
    set(groot, 'DefaultLegendInterpreter', 'latex');

    %% Load in color map
    mapObj = load("fine_red_blue_cmap.mat");
    cmap = mapObj.cmap;
    blueCol = cmap(1, :);
    redCol = cmap(end, :);
    blackCol = [0, 0, 0];

    %% Parameters
    tmax = 1;
    ts = linspace(1e-6, tmax, 1e3)';
    freq = 100;

    % Spatial limits
    xMin = 0;
    xMax = 3;
    zMax = 0.8 + 1e-10;

    %% Substrate functions
    if dimension == "2D"
        types = ["stationary", "flat", "curved"];
        colors = [blackCol; redCol; blueCol];
        displayNames = ["Free-surface: Stationary substrate", ...
            "Free-surface: Flat substrate", "Free-surface: Curved substrate"];
        dirName = "Two-dimensional_Figures";
    elseif dimension == "axi"
        types = ["stationary", "flat"];
        colors = [blackCol; redCol];
        displayNames = ["Free-surface: Stationary substrate", ...
            "Free-surface: Flat substrate"];
        dirName = "Axisymmetric_Figures";
    else
        error("Invalid dimension");
    end

    %% Loop over time
    tiledlayout(length(types), 1);

    tileNo = 1;

    for typeIdx = 1 : length(types)

        % Set up tile
        nexttile(tileNo);
        hold on;
        tileNo = tileNo + 1;

        % Loop over time
        for tIdx = 1 : freq : length(ts)
            t = ts(tIdx);

            % Load in substrate functions
            SubstrateFunctions = substratefunctions(types(typeIdx), dimension);

            % Find minimum tau
            d = SubstrateFunctions.d;
            d_t = SubstrateFunctions.d_t;
            zeroFun = @(tau) xMax - 2 * d_t(tau) * (t - tau) - d(tau);
            tauMin = fsolve(zeroFun, 1e-6);

            % Find taus
            taus = linspace(tauMin, t, 1e2);

            % Find free-surface
            [xBars, hBars] = freesurface(t, taus, SubstrateFunctions);

            % Plot free-surface
            h(1) = plot(xBars, hBars, 'linewidth', 2, 'color', colors(typeIdx, :), ...
                'Displayname', displayNames(typeIdx));
            hold on;

        end

        % Plot turnover evolution
        tLongs = linspace(0, 2 * tmax);
        ds = SubstrateFunctions.d(tLongs);
        Js = SubstrateFunctions.J(tLongs);
        plot(ds, Js, 'color', colors(typeIdx, :), 'linestyle', '--');

        %% Figure settings

        % Axes labels
        if dimension == "2D"
            xlabel('$\bar{x}$');
            height = 600;
        else
            xlabel('$\bar{r}$');
            height = 400;
        end
        ylabel('$\bar{z}$');

        box on;
        grid on;

        % Axes limits
        xlim([xMin, xMax]);
        ylim([0, zMax]);

        % Legend
        legend(h(1), 'location', 'northwest');

    end

    set(gcf,'position', [0, 0, 800, height]);

    %% Create figures
    filename = sprintf("JetEvolution_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);

    %% Function definitions
    function [xBars, hBars] = freesurface(t, taus, SubstrateFunctions)

        ds = SubstrateFunctions.d(taus);
        d_ts = SubstrateFunctions.d_t(taus);
        d_tts = SubstrateFunctions.d_tt(taus);
        Js = SubstrateFunctions.J(taus);

        xBars = 2 * d_ts .* (t - taus) + ds;
        hBars = (d_ts .* Js) ./ (d_ts - 2 * d_tts .* (t - taus));
    end

end