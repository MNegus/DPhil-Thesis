%% cSolutionComparison
% Compares the exact solution to c(t) by solving
% exp(2 * c) + 4 * exp(c) + 2 * c + 1 = pi * d / (epsilon^2 * J),
% to the asymptotic solution
% c = 0.5 * log(pi * d / (epsilon^2 * J))

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

% Times 
ts = linspace(1e-3, 10, 1e3);

%% Solutions for c
% Asymptotic solution
c_leading = 2 * sqrt(pi * d(ts) ./ (4 * epsilon^2 * J(ts)));
c_first = 2 * sqrt(1 + pi * d(ts) ./ (4 * epsilon^2 * J(ts))) - 2;

% Numerical solution
zero_fun = @(c) c.^2 + 4 * c + 2 * log(c) + 1 ...
    - pi * d(ts) ./ (epsilon^2 * J(ts));
c_num = fsolve(zero_fun, c_leading);

% Plot and compare
figure(1);
hold on;
plot(ts, c_num, 'Linestyle', '--', 'linewidth', 2);
plot(ts, c_leading, 'linewidth', 2);
plot(ts, c_first, 'linewidth', 2);
set(gca, 'Yscale', 'log');
legend('Numerical', 'Leading', 'First');
%% Solutions for inner force
% Finner = @(c) 8 * epsilon * d_t(ts) .* J(ts) .* c / pi;
FinnerFull = (16 * d_t(ts).^2 / pi) ...
    .* (sqrt(pi * d(ts) .* J(ts) / 4 + epsilon^2 * J(ts).^2) - epsilon * J(ts));
FinnerFirst = (16 * d_t(ts).^2 / pi) ...
    .* (sqrt(pi * d(ts) .* J(ts) / 4));
FinnerSecond = (16 * d_t(ts).^2 / pi) ...
    .* (sqrt(pi * d(ts) .* J(ts) / 4)  - epsilon * J(ts) );


figure(2);
hold on;
plot(ts, FinnerFirst);
plot(ts, FinnerSecond);
plot(ts, FinnerFull);
legend('First', 'Second', 'Full');

