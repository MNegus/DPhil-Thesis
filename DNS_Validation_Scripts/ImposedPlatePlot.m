%% ImposedPlatePlot.m
% Script to plot the time dependence of the imposed plate motion used for
% the validation.

clear;
close all;

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Parameters
% Time parameters used by the simulation
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
zeroIdx = sum(ts < 0);

a = 0.00125; % Imposed coefficient for the plate motion
omega = 12; % Frequency of the cosine

%% Plate solution
w = @(t) a * (1 - cos(omega *  t) + t.^2); % Function for w for t > 0

ws = zeros(size(ts));
ws(1 : zeroIdx) = 0;
size(ws(zeroIdx + 1 : end))
size(w(ts_analytical))

ws(zeroIdx + 1 : end) = w(ts_analytical);

%% Plot the solution
figure(1);
plot(ts, ws, 'linewidth', 2, 'color', 'black');

%% Figure options
% Axes labels
xlabel("$t$");
ylabel("$w(t)$");

xlim([-IMPACT_TIME, 0.7]);
grid on;

% Set figure position
x0=400;
y0=400;
height=300;
width=800;
set(gcf,'position',[x0,y0,width,height]);
set(gcf, 'Renderer', 'Painters');

% Save figure
exportgraphics(gcf, "dns_validation_figures/ImposedPlate.png", "Resolution", 300);
savefig(gcf, "dns_validation_figures/ImposedPlate.fig");
