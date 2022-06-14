%% PressureForceComparison
% Plots both the pressure and force along the plate, comparing the outer,
% inner and composite solutions

clear;
close all;

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Parameter definitions
epsilon = 0.01;
% Time dependent quantities
d = @(t) 2 * sqrt(t);
d_t = @(t) 1 ./ sqrt(t);
J = @(t) pi * d(t) ./ (8 * d_t(t).^2);
A = @(t) d(t) .* d_t(t);
C = @(t) A(t);

% Times 
ts = linspace(1e-3, 1, 1e3);

%% Force values
% Outer force
Fouter = pi * A(ts);

% Inner force, both the full and leading-order
FinnerFull = (16 * d_t(ts).^2 / pi) ...
    .* (sqrt(pi * d(ts) .* J(ts) / 4 + epsilon^2 * J(ts).^2) - epsilon * J(ts));
FinnerLeading = (16 * d_t(ts).^2 / pi) ...
    .* (sqrt(pi * d(ts) .* J(ts) / 4));

% Overlap force
Foverlap = 2 * sqrt(2) * C(ts);

% Composite
FcompFull = Fouter + FinnerFull - Foverlap;
FcompLeading = Fouter + FinnerLeading - Foverlap;

Foverlap - FinnerLeading

%% Compare forces
figure(1);
hold on;
% plot(ts, Fouter, 'color', 'black', 'linewidth', 2, 'linestyle', '--');
% plot(ts, FinnerLeading, 'color', 'black', 'linewidth', 2, 'linestyle', ':');  
% plot(ts, FinnerFull, 'color', 0.5 * [1 1 1], 'linewidth', 2, 'linestyle', ':');
% plot(ts, Foverlap, 'color', 'red');
plot(ts, FcompLeading, 'color', 'black', 'linewidth', 2);
plot(ts, FcompFull, 'color', 0.5 * [1 1 1], 'linewidth', 2); 

% legend('Outer', 'Inner: Leading', 'Inner: Full', 'Composite: Leading', 'Composite: Full');
legend('Composite: Leading', 'Composite: Full');