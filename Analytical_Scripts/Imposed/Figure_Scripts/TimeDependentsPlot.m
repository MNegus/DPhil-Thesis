function TimeDependentsPlot(dimension)
%TIMEDEPENDENTSPLOT Function to plot the time-dependent quantities
%   Script to make figures comparing the time-dependent quantities (i.e.
%   substrate position, turnover point, jet thickness etc. for either 2D
%   (where dimension == '2D') or axisymmetric (where dimension == 'axi').

    close all;

    addpath("../");
    addpath("../Forces");
    addpath("../Energies");

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

    %% Substrate functions
    if dimension == "2D"
        types = ["stationary", "flat", "curved"];
        colors = [blackCol; redCol; blueCol];
        displayNames = ["Stationary substrate", "Flat substrate", "Curved substrate"];
        dirName = "Two-dimensional_Figures";
    elseif dimension == "axi"
        types = ["stationary", "flat"];
        colors = [blackCol; redCol];
        displayNames = ["Stationary substrate", "Flat substrate"];
        dirName = "Axisymmetric_Figures";
    else
        error("Invalid dimension");
    end
    
    % Load substrate functions
    for typeIdx = 1 : length(types)
        SubstrateFunctions(typeIdx) = substratefunctions(types(typeIdx), dimension);
    end

    %% Load parameters
    % Load epsilon
    [epsilon, ~, ~, ~] = substrateparameters(); 

    % Set times
    tmax = 1;
    ts = linspace(0, tmax, 1e3)';

    %% Substrate motion plot
    if dimension == "2D"
        fig = tiledlayout(1, 3);

        % Need the flat functions struct
        FlatFunctions = substratefunctions('flat', dimension);

        % Substrate position
        nexttile;
        ws = FlatFunctions.a(ts);
        plot(ts, ws, 'color', redCol, 'linewidth', 2);
        xlabel("$t$");
        ylabel("$w(t$)");
        grid on;

        % Substrate velocity
        nexttile;
        w_ts = FlatFunctions.a_t(ts);
        plot(ts, w_ts, 'color', redCol, 'linewidth', 2);
        xlabel("$t$");
        ylabel("$\dot{w}(t)$");
        grid on;

        % Substrate acceleration
        nexttile;
        w_tts = FlatFunctions.a_tt(ts);
        plot(ts, w_tts, 'color', redCol, 'linewidth', 2);
        xlabel("$t$");
        ylabel("$\ddot{w}(t)$");
        grid on;

        % Figure settings
        set(gcf,'position', [100, 100, 800, 250]);
        pause(0.1);

        % Export figure
        filename = sprintf("FlatImposed_%s", dimension);
        savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
        exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
        exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);
    end

    %% Curved substrate spatial evolution
    figNo = 2;
    if dimension == "2D"
        figure(figNo);
        figNo = figNo + 1;
        hold on;
        CurvedFunctions = substratefunctions('curved', dimension);
        
        L = epsilon * 2 * sqrt(tmax);
        xs = linspace(0, L, 1e3);
        
        % Determine number of times
        freq = 100;
        tIdxs = 1 : freq : length(ts);
        
        length(cmap)
        length(tIdxs)
        colorFreq = floor(length(cmap) / length(tIdxs));
        
        for k = 1 : length(tIdxs)
            t = ts(tIdxs(k));
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
    
        % Time arrow
        annotation('arrow',[0.242 0.242],...
        [0.709 0.873]);

        % Time textbox
        annotation('textbox',...
            [0.250180327868853 0.84 0.0826065573770493 0.0700000000000006],...
            'String','$t$',...
            'LineStyle','none',...
            'Interpreter','latex',...
            'FitBoxToText','off', ...
            'Fontsize', 18);

        % Export figure
        filename = sprintf("CurvedImposed_%s", dimension);
        savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
        exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
        exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);
    end
    

    
    %% Turnover point plot
    figure(figNo);
    figNo = figNo + 1;
    hold on;
    for typeIdx = 1 : length(types)
        plot(ts, SubstrateFunctions(typeIdx).d(ts), 'color', colors(typeIdx, :), ...
            'linewidth', 2);
    end
    legend(displayNames, 'Location', 'southeast');
    xlabel("$t$");
    ylabel("$d_0(t)$");
    grid on;

    % Figure settings
    set(gcf,'position', [100, 100, 800, 400]);
    pause(0.1);

    % Export figure
    filename = sprintf("TurnoverPoints_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);
    
    %% Turnover derivative plot
    figure(figNo);
    figNo = figNo + 1;
    hold on;
    for typeIdx = 1 : length(types)
        plot(ts, SubstrateFunctions(typeIdx).d_t(ts), 'color', colors(typeIdx, :), ...
            'linewidth', 2);
    end
    legend(displayNames, 'Location', 'northeast');
    xlabel("$t$");
    ylabel("$\dot{d}_0(t)$");
    grid on;

    ylim([0, 10]);
    
    % Figure settings
    set(gcf,'position', [100, 100, 800, 400]);
    pause(0.1);

    % Export figure
    filename = sprintf("TurnoverPointDerivativess_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);

    %% Jet thickness
    figure(figNo);
    figNo = figNo + 1;
    hold on;
    for typeIdx = 1 : length(types)
        plot(ts, SubstrateFunctions(typeIdx).J(ts), 'color', colors(typeIdx, :), ...
            'linewidth', 2);
    end

    legend(displayNames, 'Location', 'northwest');
    xlabel("$t$");
    ylabel("$J(t)$");
    grid on;

    % Figure settings
    set(gcf,'position', [100, 100, 800, 400]);
    pause(0.1);

    % Export figure
    filename = sprintf("JetThickness_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);

    %% Force plot

    figure(figNo);
    figNo = figNo + 1;
    hold on;

    for typeIdx = 1 : length(types)
        % Load forces
        [Fs_composite, Fs_outer, ~] ...
            = substrateforce(ts, SubstrateFunctions(typeIdx));
        
        % Plot composite
        h(typeIdx) = plot(ts, Fs_composite, 'color', colors(typeIdx, :), ...
            'linewidth', 2, 'Displayname', displayNames(typeIdx));
       
        % Plot outer
        plot(ts, Fs_outer, 'color', colors(typeIdx, :), ...
            'linewidth', 2, 'linestyle', '--');
    end

    legend(h(1 : length(types)), 'Location', 'northwest');

    xlabel("$t$");
    ylabel("$F(t)$");

    if dimension == "2D"
        ylim([0, 15]);
    else
        ylim([0, 2]);
    end
    grid on;

    % Figure settings
    set(gcf,'position', [100, 100, 800, 400]);
    pause(0.1);

    % Export figure
    filename = sprintf("Force_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);

    %% Energy plot
    figure(figNo);
    figNo = figNo + 1;
    hold on;

    for typeIdx = 1 : length(types)
        % Load energies
        [Es_outer, Es_jets] = dropletenergy(ts, SubstrateFunctions(typeIdx));
        % Plot outer
        h(typeIdx) = plot(ts, Es_outer, 'color', colors(typeIdx, :), ...
            'linewidth', 1.5, 'Displayname', displayNames(typeIdx));
       
        % Plot jets
        plot(ts, Es_jets, 'color', colors(typeIdx, :), ...
            'linewidth', 3, 'linestyle', '--');
    end

    legend(h(1 : length(types)), 'Location', 'northwest');
    xlabel("$t$");
    ylabel("$E_{K}(t)$");
    grid on;
    
    if dimension == "2D"
        ylim([0, 0.06]);
    else
        ylim([0, 0.004]);
    end

    % Figure settings
    set(gcf,'position', [100, 100, 800, 400]);
    pause(0.1);

    % Export figure
    filename = sprintf("Energy_%s", dimension);
    savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
    exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
    exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);
end
