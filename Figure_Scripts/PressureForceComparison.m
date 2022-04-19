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

ts = linspace(1e-3, 10, 20);

%% Forces values
FouterExact = pi * A(ts);
FcompExact = pi * A(ts) + 8 * sqrt(J(ts) / pi) .* d(ts).^(5/2) - 2 * sqrt(2) * A(ts);

FouterNumerical = zeros(size(ts));
FcompNumerical = zeros(size(ts));

%% Plots in time
figure(1);
for k = 1 : length(ts)
    t = ts(k)

    %% Inner pressure
    c0 = log(sqrt(pi * d(t) ./ (epsilon^2 * J(t))));
    sigmas = linspace(-c0, 1e3, 1e5);
    xTildes = (J(t) / pi) * (2 * sigmas - exp(-2 * sigmas) - 4 * exp(-sigmas) - 1);
    xs = epsilon * d(t) + epsilon^3 * xTildes;

    psInner = (1 / epsilon^2) * 2 * d_t(t)^2 * exp(-sigmas) ./ (1 + exp(-sigmas)).^2;

    %% Outer and overlap pressure
    xHats = xs / epsilon;
    xHats = xHats(xHats < d(t));
    size(xHats)
    psOuter = A(t) ./ (epsilon * sqrt(d(t).^2 - xHats.^2));
    psOverlap = C(t) ./ (epsilon * sqrt(2 * d(t)) * sqrt(d(t) - xHats));
    
    %% Composite pressure
    psComp = psInner;
    psComp(xHats < d(t)) = psComp(xHats < d(t)) + psOuter - psOverlap;

    %% Numerical calculation
    FouterNumerical(k) = 2 * trapz(epsilon * xHats, psOuter);
    FcompNumerical(k) = 2 * trapz(xs, psComp);
    
    %% Plots
    plot(epsilon * xHats, psOuter);
    hold on;
    plot(xs, psInner);
    plot(xs, psComp);
    hold off;
    drawnow;
    ylim([0, 1.5 * max(psComp)]);
    pause(0.1);

end

%% Plot comparisons
close(figure(2));
figure(2);
hold on;
plot(ts, FouterNumerical, 'Color', 'Black', 'Linewidth', 1);
plot(ts, FouterExact, 'Color', 'Black', 'Linewidth', 1, 'Linestyle', '--');
plot(ts, FcompNumerical, 'Color', 0.5 * [1 1 1], 'Linewidth', 1);
% plot(ts, FcompExact, 'Color', 0.5 * [1 1 1], 'Linewidth', 2, 'Linestyle', '--');

