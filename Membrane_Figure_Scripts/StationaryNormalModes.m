%% StationaryNormalModes
% Validates the normal modes solution when the membrane is stationary


clear;
close all;


addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/Pressures");
addpath("../Analytical_Scripts/MembraneSolution/FiniteDifference/");
addpath("../Analytical_Scripts/MembraneSolution/FiniteDifference/PressuresFD/");
addpath("../Analytical_Scripts/MembraneSolution/NormalModes/");
addpath("../Analytical_Scripts/ImposedSolution/")

parent_dir = "/home/michael/scratch/AnalyticalMembraneTests/";

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;

redCol = cmap(end, :);
blueCol = cmap(1, :);

%% Figure options
fontsize = 10;
lineWidth = 1.25;
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize', fontsize);
set(0, 'DefaultTextFontSize', fontsize);
set(0,'defaultLegendFontSize', fontsize, 'DefaultLegendFontSizeMode','manual');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(0, 'DefaultLegendInterpreter', 'latex');

% Scatter size
sz = 20;

% Tiled layout sizes
width = 6;
height = 5;

%% Parameters
EPSILON = 0.1;
L = 16;
T_MAX = 0.25;
DELTA_T = 1e-4;

% FD parameters
N_MEMBRANE = 21848;

DELTA_X = L / (N_MEMBRANE - 1); 
M = N_MEMBRANE - 1;
xs = (0 : DELTA_X : L - DELTA_X)';

ts_analytical = 0 : DELTA_T : T_MAX;

tTest = 1 / 16;
tTestIdx = tTest / DELTA_T;

%% Pressure solution

% Exact solutions
d = 2 * sqrt(tTest);
d_t = 1 ./ sqrt(tTest);

xHats = linspace(0, d, 1e4 + 1);
xHats = xHats(1 : end - 1);

pExact = 2 ./ sqrt(d^2 - xHats.^2);


% Normal modes range
% N_Ms = 10.^(1 : 6);
N_Ms = 2.^(1 : 17);

% Errors
errors = zeros(size(N_Ms));

% Plot pressure
tiledlayout(1, 3);
height = 2.5;
width = 6;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

nexttile([1,2])
hold on;

LineColorIdxs = floor(linspace(1, length(cmap), 5));
lineCount = 1;

for NIdx = 1 : length(N_Ms)
    N_M = N_Ms(NIdx)
    
    % Define ks
    ks = pi * (2 * (1 : N_M) - 1) / (2 * L);

    % Determine ps
    ps = (2 * pi * EPSILON / L) * besselj(0, EPSILON * ks * d) * cos(EPSILON * ks' * xHats);

    % Determine error
    errors(NIdx) = sqrt(sum((ps - pExact).^2 ./ pExact.^2) / length(pExact)) ;

    % Plot pressure
    if (mod(NIdx - 1, 4) == 0)
        displayName = append("$N_M = ", num2str(N_M), "$");
        lineColor = cmap(LineColorIdxs(lineCount), :);
        plot(xHats, ps, 'DisplayName', displayName, 'Color', lineColor, ...
            'LineWidth', lineWidth);
        lineCount = lineCount + 1;
    end

end

% Plot exact
plot(xHats, pExact, 'Color', 'Black', 'LineStyle', '--', ...
    'LineWidth', lineWidth, 'DisplayName', 'Exact');

legend('Location', 'Northwest', 'NumColumns', 2);
set(gca, 'Yscale', 'log');
set(gca,'xminorgrid','off','yminorgrid','off');
ylim('padded');
xlim('padded');
xlabel("$\hat{x}$")
ylabel("$\hat{p}_0(\hat{x}, 0, t)$")
grid on;
box on;


% Plot error
nexttile(3)
hold on;
scatter(N_Ms, errors, 'Color', blueCol);

N_MsExtended = linspace(N_Ms(1), 20 * N_Ms(end));
plot(N_MsExtended, 2 * errors(end) * sqrt(N_Ms(end) ./ N_MsExtended), ...
    'LineStyle', ':', 'Color', 'black', 'LineWidth', lineWidth)

set(gca, 'Yscale', 'log');
set(gca, 'XScale', 'log');
ylim([10^-2, 10^0.5]);
xlim([0.5, 10^6.5]);
xticks([10^0, 10^2, 10^4, 10^6]);
set(gca,'xminorgrid','off','yminorgrid','off');
xlabel("$N_M$")
ylabel("Error")
grid on;
box on;

% Pause
pause(0.01);

% Output
figname = "MembraneFigures/StationaryNMPressure";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);

%% Energy distribution
tTest = 1;
ts = linspace(0, tTest, 1e3);
ds = 2 * sqrt(ts);
EExact = EPSILON^2 * pi * ts;

% Normal modes
N_Ms = 2.^(1 : 13);

% Errors
errors = zeros(size(N_Ms));

% Plot energy
tiledlayout(1, 3);
height = 3;
width = 6;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

LineColorIdxs = floor(linspace(1, length(cmap), 5));
lineCount = 1;

nexttile([1, 2]);
hold on;

for NIdx = 1 : length(N_Ms)
    N_M = N_Ms(NIdx);

    Es = zeros(size(ts));

    for n = 1 : N_M
        k = pi * (2 * n - 1) / (2 * L);
        Es = Es + besselj(1, EPSILON * k * ds) .* sin(EPSILON * k * ds) / k^2;
    end

    Es = (EPSILON * pi * ds / L) .* Es;

    % Determine error
    errors(NIdx) = sqrt(sum((Es(2 : end) - EExact(2 : end)).^2 ./ EExact(2 : end).^2) / length(EExact));


    if (mod(NIdx - 1, 3) == 0)
        displayName = append("$N_M = ", num2str(N_M), "$");
        lineColor = cmap(LineColorIdxs(lineCount), :);
        plot(ts, Es, 'DisplayName', displayName, 'Color', lineColor, ...
            'LineWidth', lineWidth);
        lineCount = lineCount + 1;
    end
    

    drawnow;
    pause(0.01);
end

plot(ts, EExact, 'Color', 'Black', 'LineStyle', '--', ...
    'LineWidth', lineWidth, 'DisplayName', 'Exact');
legend('Location', 'Northoutside', 'NumColumns', 2);
ylim('padded');
xlim('padded');
xlabel("$t$")
ylabel("$E_{K, outer}(t)$")
grid on;
box on;


% Plot error
nexttile(3)
hold on;
scatter(N_Ms, errors, 'Color', blueCol);

N_MsExtended = linspace(N_Ms(1), 20 * N_Ms(end));
plot(N_MsExtended, 2 * errors(end) * (N_Ms(end) ./ N_MsExtended).^1.5, ...
    'LineStyle', ':', 'Color', 'black', 'LineWidth', lineWidth)


set(gca, 'Yscale', 'log');
set(gca, 'XScale', 'log');
ylim([10^-4, 10^0.5]);
xlim([0.5, 10^4.5]);
xticks([10^0, 10^2, 10^4, 10^6]);
set(gca,'xminorgrid','off','yminorgrid','off');
xlabel("$N_M$")
ylabel("Error")
grid on;
box on;

% Output
figname = "MembraneFigures/StationaryNMEnergy";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);