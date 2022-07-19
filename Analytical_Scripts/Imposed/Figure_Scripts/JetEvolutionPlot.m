function JetEvolutionPlot(dimension)
%JETEVOLUTIONPLOT
% Script to plot the time evolution of the jets/splash sheet

    close all;

    addpath("../");

    %% Figure options
    set(0,'defaultTextInterpreter','latex'); %trying to set the default
    set(0,'defaultAxesFontSize', 18);
    set(0,'defaultLegendFontSize', 18, 'DefaultLegendFontSizeMode','manual');
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
        displayNames = ["Stationary substrate", ...
            "Moving substrate (rigid)", "Moving substrate (curved)"];
        dirName = "Two-dimensional_Figures";
    elseif dimension == "axi"
        types = ["stationary", "flat"];
        colors = [blackCol; redCol];
        displayNames = ["Stationary substrate", ...
            "Moving substrate"];
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
            d_tt = SubstrateFunctions.d_tt;
            zeroFun = @(tau) xMax - 2 * d_t(tau) * (t - tau) - d(tau);
            tauMin = fsolve(zeroFun, 1e-6);

            % Find taus
            taus = linspace(tauMin, t, 1e2);

            % Determinant of Jacobian, throw error if zero at any point
            detJ = d_t(taus) - 2 * d_tt(taus) .* (t - taus);
            if (sum(detJ <=0) > 1)
                detJ
                error("Singular Jacobian");
            end
            
            % Find free-surface
            [xBars, hBars, uBars] = freesurface(t, taus, SubstrateFunctions);
            
            % Find pressure (DOESNT WORK FOR AXI)
%             zBar = 0;
%             pBars = pressure(t, xBars, zBar, uBars, hBars, SubstrateFunctions);

            % Plot free-surface
            h(1) = plot(xBars, hBars, 'linewidth', 2, 'color', colors(typeIdx, :), ...
                'Displayname', displayNames(typeIdx));
%             h(1) = plot(xBars, pBars, 'linewidth', 2, 'color', colors(typeIdx, :), ...
%                 'Displayname', displayNames(typeIdx));
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
            ylabel('$\bar{h}_0(\bar{x}, t)$');
            height = 600;
        else
            xlabel('$\bar{r}$');
            ylabel('$\bar{h}_0(\bar{r}, t)$');
            height = 400;
        end
        

        box on;
        grid on;

        % Axes limits
        xlim([xMin, xMax]);
        ylim([0, zMax]);

        % Legend
        legend(h(1), 'location', 'northwest');

    end

    set(gcf,'position', [0, 0, 1000, height]);

    %% Create figures
    filename = sprintf("JetEvolution_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);

    %% Function definitions
    function [xBars, hBars, uBars] = freesurface(t, taus, SubstrateFunctions)

        ds = SubstrateFunctions.d(taus);
        d_ts = SubstrateFunctions.d_t(taus);
        d_tts = SubstrateFunctions.d_tt(taus);
        Js = SubstrateFunctions.J(taus);

        xBars = 2 * d_ts .* (t - taus) + ds;
        hBars = (d_ts .* Js) ./ (d_ts - 2 * d_tts .* (t - taus));
        uBars = 2 * d_ts;
    end

    function pBars = pressure(t, xBars, zBar, uBars, hBars, SubstrateFunctions)
        
        % Load w functions
        epsilon = SubstrateFunctions.epsilon;
        xs = epsilon * xBars;
        wBar_xxs = epsilon^2 * SubstrateFunctions.w_xx(xs, t);
        wBar_xts = epsilon * SubstrateFunctions.w_xt(xs,  t);
        wBar_tts = SubstrateFunctions.w_tt(xs, t);
        
        % Determine pressure
        pBars = -(wBar_tts + 2 * wBar_xts .* uBars + wBar_xxs .* uBars.^2) .* (hBars - zBar);
        
    end

end