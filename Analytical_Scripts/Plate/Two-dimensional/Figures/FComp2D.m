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
epsilon = 1;
% Time dependent quantities
d = @(t) 2 * sqrt(t);
d_t = @(t) 1 ./ sqrt(t);
J = @(t) pi * d(t) ./ (8 * d_t(t).^2);
A = @(t) d(t) .* d_t(t);
C = @(t) A(t);

% Times 
tmax = 0.8 - 0.125;
ts = linspace(1e-3, tmax, 1e3);

%% Solution for c
% Asymptotic solution
c_leading = 2 * sqrt(pi * d(ts) ./ (4 * epsilon^2 * J(ts)));

% Numerical solution
zero_fun = @(c) c.^2 + 4 * c + 2 * log(c) + 1 ...
    - pi * d(ts) ./ (epsilon^2 * J(ts));
c = fsolve(zero_fun, c_leading);

%% Force values
% Outer force
Fouter = pi * A(ts);

% Inner force, both the full and leading-order
Finner = 8 * epsilon * d_t(ts).^2 .* J(ts) .* c / pi;

% Overlap force
Foverlap = 2 * sqrt(2) * C(ts);

% Composite
Fcomp = Fouter + Finner - Foverlap;

%% Numerical solution
outputMat = readmatrix("output.txt");
tsNum = outputMat(:, 1) - 0.125;
Fnum = outputMat(:, 2);

%% Compare forces
figure(1);
hold on;
plot(ts, Fouter, 'color', 'black', 'linewidth', 2, 'linestyle', '--');
% plot(ts, Finner, 'color', 'black', 'linewidth', 2, 'linestyle', ':');  
% plot(ts, FinnerFull, 'color', 0.5 * [1 1 1], 'linewidth', 2, 'linestyle', ':');
% plot(ts, Foverlap, 'color', 'red');
plot(ts, Fcomp, 'color', 'black', 'linewidth', 2);
scatter(tsNum, Fnum);

% legend('Outer', 'Inner: Leading', 'Inner: Full', 'Composite: Leading', 'Composite: Full');
legend('Outer', 'Comp', 'DNS');