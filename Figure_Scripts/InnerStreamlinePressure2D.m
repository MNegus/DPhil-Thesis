%% OuterStreamlinePressure
% Script to produce a figure of the streamlines and pressure in theinner
% region for the 2D impact case. Appears in Section 3.3.4, under the Wagner
% theory chapter. 
% 
% For the figure, we neglect any motion of the membrane, taking w == 0. 
%
% ADD DESCRIPTION AS TO WHAT WE ARE DOING IN TERMS OF PLOTTING STREAMLINES
% ETC

clear;
close all;

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Load in color map
mapObj = load("red_blue_cmap.mat");
cmap = mapObj.cmap;

%% Parameter definitions
% Time dependent quantities
d = 1; % Let the turnover point be at 1
t = d^2 / 4; % Solves for t as a function of d
d_t = 1 / sqrt(t); % Time derivative of d
J = pi * d / (8 * d_t^2); % Jet thickness

% Anonymous functions
tzeta = @(eta) (J / pi) * (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);

%% Figure properties
figure1 = figure();
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Set colour map to cool-warm
colormap(cmap);

% Axes labels
xlabel('$\hat{x}$');
ylabel('$\hat{z}$');

% Show axes grid
grid on;


%% Plot the free-surface
noPoints = 1e3;
xiMax = 100;

% Lower values
xisLower = linspace(-xiMax, 0, noPoints / 2); 
xsLower = (J / pi) * (exp(xisLower) - xisLower - 1);
hsLower = (J / pi) * (pi + 4 * exp(xisLower / 2));

% Upper values
xisUpper = linspace(0, xiMax, noPoints / 2); 
xsUpper = (J / pi) * (xisUpper - log(1 + xisUpper));
hsUpper = (J / pi) * (pi + 4 * sqrt(xisUpper + 1));

% Combine and plot
plot([xsLower, xsUpper], [hsLower, hsUpper], 'color', 'black', 'linewidth', 2);

%% Plot a selection of streamlines
% thetas = linspace(-pi, pi, 1e3);


% kappas = linspace(0, 1.5 * pi, 10);
kappas = linspace(0, pi, 10);

for kappa = kappas(1:10)
    
    % Solve for thetaSwitch, the point at which r < 1
    thetaSwitch = fsolve(@(theta) kappa - theta - sin(theta), 0.5 * pi);
    
    
    % theta from 0 to kappa
    thetasUpper = linspace(0, thetaSwitch, 1e3);
%     thetasLower = linspace(thetaSwitch, kappa, 1e3);
    idxs = 1 ./ 2.^linspace(0, 1000, 1e3);
    thetasLower = thetaSwitch * idxs + kappa *  (1 - idxs);
    thetas = [thetasUpper, thetasLower];

    rs = (kappa - thetas) ./ sin(thetas); % Finds solution for rs

    % Creates array of eta values from the non-negative r values
    etas = rs(rs >= 0) .* exp(1i * thetas(rs >= 0));

    % Plots etas
    % figure(2);
    % plot(real(etas), imag(etas));

    % Find tzeta values
    tzetas = tzeta(etas);

    % Plot the streamline
    figure(1);
    plot(real(tzetas), imag(tzetas));
    
    xlim([-2, 2]);
    ylim([0, 2]);
    drawnow;
    
    figure(2);
    plot(real(etas), imag(etas));
    hold on;
    xlim([-2, 2]);
    ylim([0, 2]);
    pause(0.1);
end



%% Figure scalings
xlim([-2, 2]);
ylim([0, 2]);



