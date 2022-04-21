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

ts = linspace(1e-3, 50, 50);

%% Solve for c0 function
c0Guess = 0.5 * log(pi * d(ts) ./ (epsilon^2 * J(ts)));
zero_fun = @(c) exp(2 * c) + 4 * exp(c) + 2 * c + 1 - pi * d(ts) ./ (epsilon^2 * J(ts));
c0s = fsolve(zero_fun, c0Guess);
% c0s = c0Guess

%% Forces values
FouterExact = pi * A(ts);
% FinnerExact = 8 * sqrt(J(ts) / pi) .* d(ts).^(5/2
% FinnerExact = 8 * epsilon * d_t(ts).^2 .* J(ts) .* exp(c0s) / pi;
FinnerExact = 8 * d_t(ts).^2 .* sqrt(d(ts) .* J(ts)) / pi;
FcompExact = FouterExact + FinnerExact - 2 * sqrt(2) * A(ts);

FouterNumerical = zeros(size(ts));
FinnerNumerical = zeros(size(ts));
FcompNumerical = zeros(size(ts));

%% Plots in time
figure(1);
for k = 1 : length(ts)
    t = ts(k)

    %% Inner pressure
    
    sigmas = linspace(-c0s(k), 1e3, 1e5);
    xTildes = (J(t) / pi) * (2 * sigmas - exp(-2 * sigmas) - 4 * exp(-sigmas) - 1);
    xs = epsilon * d(t) + epsilon^3 * xTildes;

    psInner = (1 / epsilon^2) * 2 * d_t(t)^2 * exp(-sigmas) ./ (1 + exp(-sigmas)).^2;

    %% Outer and overlap pressure
    xHats = xs / epsilon;
    xHats = xHats(xHats < d(t));
    psOuter = A(t) ./ (epsilon * sqrt(d(t).^2 - xHats.^2));
    psOverlap = C(t) ./ (epsilon * sqrt(2 * d(t)) * sqrt(d(t) - xHats));
    
    %% Composite pressure
    psComp = psInner;
    psComp(xHats < d(t)) = psComp(xHats < d(t)) + psOuter - psOverlap;

    %% Numerical calculation
    
    FouterNumerical(k) = 2 * trapz(xHats, epsilon * psOuter);
    FinnerNumerical(k) = 2 * trapz(xs, psInner);
    
    FcompNumerical(k) = 2 * trapz(xs, psComp);
    
    %% Plots
    plot(epsilon * xHats, psOuter);
    hold on;
    plot(xs, psInner);
    plot(xs, psComp, 'Color', 'Black', 'linewidth', 1);
    hold off;
    
    xlim([-0.01, 0.2]);
    ylim([0, 700]);
    drawnow;
    pause(0.1);

end

%% Plot comparisons
close(figure(2));
figure(2);
hold on;
plot(ts, FouterNumerical, 'Color', 'Black', 'Linewidth', 2, 'Displayname', 'Numerical Outer');
plot(ts, FouterExact, 'Color', 'Black', 'Linewidth', 2, 'Linestyle', '--', 'Displayname', 'Exact Outer');

plot(ts, FinnerNumerical, 'Color', 'Blue', 'Linewidth', 2, 'Displayname', 'Numerical Inner');
plot(ts, FinnerExact, 'Color', 'Blue', 'Linewidth', 2, 'Linestyle', '--', 'Displayname', 'Exact Inner');

plot(ts, FcompNumerical, 'Color', 0.5 * [1 1 1], 'Linewidth', 1);
plot(ts, FcompExact, 'Color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--', 'Displayname', 'Composite');
legend();
% ylim([5, 10]);

