%% OuterInnerSurfaceComposite.m
% Plots the composite solution between the outer and inner region

clear;
close all;

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

figure(1);
hold on;

% Axes labels
xlabel('$x$');
ylabel('$z$');

box on;
grid on;

%% Parameters
epsilon = 0.01;

d = 0.1; % Turnover point
t = d^2 / (epsilon^2 * 4); % Time
d_t = 1 / sqrt(t); % Time derivative of d
J = pi * d / (8 * d_t^2); % Jet thickness

xMax = 0.5; % Maximum value of x

%% Anonymous functions
x = @(xi) epsilon * d + (epsilon^3 * J / pi) * (xi - log(1 + xi));

%% Inner solution
% Solve for xiMax
xiGuess = (pi / (epsilon^3 * J)) * (xMax - epsilon * d);
zeroFun = @(xi) x(xi) - xMax;
xiMax = fzero(zeroFun, xiGuess);

% Create xis
xis = linspace(0, xiMax, 1e4);

% Inner solution
xs = x(xis);
hsInner = (J / pi) * (pi + 4 * sqrt(xis + 1));

%% Outer solution 
xHats = xs / epsilon;
hsOuter = 0.5 * xHats .* sqrt(xHats.^2 - d^2);

%% Overlap solution
hsOverlap = d^(3/2) * sqrt(xHats - d) / sqrt(2);

%% Composite solution
hsComposite = epsilon^2 * hsOuter + epsilon^3 * hsInner - epsilon^2 * hsOverlap;
hsOuter(1)
hsInner(1)
hsOverlap(1)
hsComposite(1)

%% Plot free-surface solutions
plot(xs, epsilon^2 * hsOuter);
plot(xs, epsilon^3 * hsInner);
plot(xs, epsilon^2 * hsOverlap);
plot(xs, hsComposite);
% plot(xs, hsUndisturbed, 'linestyle', '--', 'Color', 'black');

legend(["Outer", "Inner", "Overlap", "Composite"]);


% plot(xis, x(xis) - xMax)


