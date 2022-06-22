%% JetStreamlinePressure2D.m
% Script to produce a figure of the streamlines and pressure in the jet
% region for the 2D impact case. Appears in Section 3.3.5, under the Wagner
% theory chapter. 
% 
% For the figure, we neglect any motion of the membrane, taking w == 0. 
%

clear;
close all;

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
[epsilon, L, q, omega] = quadraticparameters(); 

tmax = 1;
ts = linspace(1e-6, tmax, 100);

%% Loop over time
figure(1);
for t = ts

    % Stationary substrate coefficients
    zeroTerm = zeros(size(t));
    StationarySubstrateCoefficients ...
        = substratecoefficients(zeroTerm, zeroTerm, zeroTerm, zeroTerm, zeroTerm, zeroTerm, epsilon);

    % Flat substrate coefficients
    [aFlats, a_tFlats, a_ttFlats, bFlats, b_tFlats, b_ttFlats] ...
        = flatsubstrate(t, q, omega); 
    FlatSubstrateCoefficients ...
        = substratecoefficients(aFlats, bFlats, a_tFlats, b_tFlats, a_ttFlats, b_ttFlats, epsilon);

    % Quadratic substrate coefficients
    [aQuads, a_tQuads, a_ttQuads, bQuads, b_tQuads, b_ttQuads] ...
        = quadraticsubstrate(t, L, q, omega);
    QuadSubstrateCoefficients ...
        = substratecoefficients(aQuads, bQuads, a_tQuads, b_tQuads, a_ttQuads, b_ttQuads, epsilon);

    %% Determine time dependents
    % Stationary substrate time dependents
    StationaryTimeDependents = timedependents(t, StationarySubstrateCoefficients);

    % Flat substrate time dependents
    FlatTimeDependents = timedependents(t, FlatSubstrateCoefficients);

    % Quadratic substrate time dependents
    QuadTimeDependents = timedependents(t, QuadSubstrateCoefficients);


    %% Spatial limits
    % Range for x
    xMin = 0;
    xMax = 5;
    zMax = epsilon;

    %Finds range for tau
    tauMin = 1e-3;
    taus = linspace(tauMin, t, 1e3)';

    %% Free-surface plotting



    % Array of types
    typeArr = ["Stationary substrate solution", ...
        "Flat substrate solution", "Quadratic substrate solution"];
    % typeArr = ["Stationary substrate solution"];


    
    for typeIdx = 1 : length(typeArr)
        type = typeArr(typeIdx);

        if type == "Stationary substrate solution"

           lineColor = 'black';

            % Load in substrate coefficients
            zeroTerm = zeros(size(taus));
            as = zeroTerm;
            bs = zeroTerm;
            SubstrateCoefficients ...
                = substratecoefficients(zeroTerm, zeroTerm, zeroTerm, ...
                    zeroTerm, zeroTerm, zeroTerm, epsilon);
        elseif type == "Flat substrate solution"

            lineColor = redCol;

            % Load in substrate coefficients
            [as, a_ts, a_tts, bs, b_ts, b_tts] ...
                = flatsubstrate(taus, q, omega);
            SubstrateCoefficients ...
                = substratecoefficients(as, bs, a_ts, b_ts, a_tts, b_tts, epsilon);
        else
            lineColor = blueCol;

            % Load in substrate coefficients
            [as, a_ts, a_tts, bs, b_ts, b_tts] ...
                = quadraticsubstrate(taus, L, q, omega);
            SubstrateCoefficients ...
                = substratecoefficients(as, bs, a_ts, b_ts, a_tts, b_tts, epsilon);
        end

        % Find free-surface
        [xs, hs] = freesurface(t, taus, SubstrateCoefficients);

        % Find substrate solution
        ws = (as(end) * ones(size(xs)) + bs(end) * xs.^2);
        
        % Plot free-surface
        plot(xs, -epsilon^2 * ws + epsilon^3 * hs, 'linewidth', 2, 'color', lineColor, ...
            'Displayname', type);
        hold on;
        
        %% Plot the substrate
        plot(xs, -epsilon^2 * ws, ...
            'Color', lineColor, 'Linestyle', ':', 'Linewidth', 2);
    end
    hold off;
    


    %% Figure properties

    % Axes labels
    xlabel('$\bar{x}$');
    ylabel('$\bar{z}$');

    box on;

%     xlim([xMin, xMax]);
    ylim([0, zMax]);
    drawnow;
    pause(0.1)
    % ylim([0, zMax]);
    % pbaspect([1 zMax / (xMax- xMin) 1]);

    % x-axis settings
    % set(gca, 'xtick',[1, 2, 3]);
    % xNames = {'$d_0(t)$'; '$2 d_0(t)$'; '$3 d_0(t)$'};
    % set(gca, 'XTickLabel', xNames);
    % 
    % % y-axis settings
    % set(gca, 'ytick',[0, 0.5 * J(t), J(t)]);
    % yNames = {'$0$'; '$0.5 \, J(t)$'; '$J(t)$'};
    % set(gca, 'YTickLabel', yNames);

    % Set figure size
    set(gcf,'position', [0, 0, 800, 400]);

end


%% Create figures
% Export png
% exportgraphics(gca,'png/JetStreamlinePressure2D.png', 'Resolution', 300);

%% Function definitions
function [xs, hs] = freesurface(t, taus, SubstrateCoefficients)

    TimeDependents = timedependents(taus, SubstrateCoefficients);
    
    ds = TimeDependents.ds;
    d_ts = TimeDependents.d_ts;
    d_tts = TimeDependents.d_tts;
    Js = TimeDependents.Js;
    
    xs = 2 * d_ts .* (t - taus) + ds;
    hs = (d_ts .* Js) ./ (d_ts - 2 * d_tts .* (t - taus));
    
    xs
    
end