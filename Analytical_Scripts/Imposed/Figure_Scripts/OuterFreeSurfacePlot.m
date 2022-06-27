function OuterFreeSurfacePlot(dimension)
    %%OuterFreeSurfacePlot
    %   Plots the free surface for the stationary, flat and quadratic substrate
    %   cases.

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
    blackCol = [0, 0, 0];

    %% Load parameters
    % Substrate parameters
    [epsilon, ~, ~, ~] = substrateparameters(); 

    tmax = 1;
    ts = [0.5, tmax];


    L = epsilon * 2 * sqrt(tmax);
    xs = linspace(0, 2 * L, 1e3);

    %% Substrate functions
%     dimension = "2D";
    if dimension == "2D"
        types = ["stationary", "flat", "curved"];
        colors = [blackCol; redCol; blueCol];
        displayNames = ["Stationary substrate", ...
            "Flat substrate", "Curved substrate"];
        dirName = "Two-dimensional_Figures";
    elseif dimension == "axi"
        types = ["stationary", "flat"];
        colors = [blackCol; redCol];
        displayNames = ["Stationary substrate", ...
            "Flat substrate"];
        dirName = "Axisymmetric_Figures";
    else
        error("Invalid dimension");
    end


    %% Free-surface plot
    tiledlayout(1, 2);

    tileNo = 1;

    
    for t = ts
        
        % Set up tile
        nexttile(tileNo);
        hold on;
        tileNo = tileNo + 1;
        
        %% Plot substrates
        for typeIdx = length(types) : -1 : 1

            % Save SubstrateCoefficients
            SubstrateFunctions = substratefunctions(types(typeIdx), dimension);
            lineColor = colors(typeIdx, :);

            if dimension == "2D"
                ws = SubstrateFunctions.w(xs, t);
            else
                ws = SubstrateFunctions.w(t) * ones(size(xs));
            end
            plot(xs / epsilon, -ws, ...
                'Color', lineColor, 'Linestyle', ':', 'Linewidth', 2);
        end
        
        %% Plot free-surfaces
        for typeIdx = length(types) : -1 : 1

            % Save SubstrateCoefficients
            SubstrateFunctions = substratefunctions(types(typeIdx), dimension);
            lineColor = colors(typeIdx, :);

            % Load turnover point
            d = SubstrateFunctions.d(t);

            % Create xs
            xHats_Free_Surface = linspace(d, 2 * L / epsilon, 1e3);

            % Determine free-surface
            hHats = outerfreesurface(xHats_Free_Surface, t, SubstrateFunctions);
            hs = hHats;

            % Plot free-surface
            h(typeIdx) = plot(xHats_Free_Surface, hs, 'linewidth', 2, 'color', lineColor, ...
                'Displayname', displayNames(typeIdx));

        end
        
        %% Figure settings
        ylim([-1, 4]);

        grid on;

        if dimension == "2D"
            xlabel("$\hat{x}$");
        else
            xlabel("$\hat{r}$");
        end
        ylabel("$\hat{z}$");
        
    end
    
    leg = legend(h(1 : length(types)), 'Orientation', 'Horizontal');
    leg.Layout.Tile = 'north';


    %% Figure settings
    % xlim([-L, L]);
    set(gcf,'position', [100, 100, 800, 400]);
    set(gcf, 'Renderer', 'painters');

    % Export figure
    filename = sprintf("OuterFreeSurface_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);

end