clear;
close all;

addpath("Pressures");
addpath("Forces");

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');


%% Parameters
epsilon = 0.1;
% tmax = 0.625 / (epsilon^2);
tmax = 1;
ts = linspace(0, tmax, 1e3)';

%% Substrate parameters, such that w(x, t) = w(t) * (L^2 - x^2);
L = 1;

% Just quadratic time dependence
% q = 0.1;
% w_coeffs = q * ts.^2;
% w_t_coeffs = 2 * q * ts;
% w_tt_coeffs = 2 * q;

% Oscillatory time dependence (with epsilons)
% q = 0.15;
% omega = 250 * epsilon^2;
% w_coeffs = q * (epsilon^4 * ts.^2 + 1 - cos(omega * ts));
% w_t_coeffs = q * (2 * epsilon^4 * ts + omega * sin(omega * ts));
% w_tt_coeffs = q * (2 * epsilon^4 + omega^2 * cos(omega * ts));

% Oscillatory time dependence (without epsilons)
q = 0.1;
omega = 2.5;
w_coeffs = q * (ts.^2 + 1 - cos(omega * ts));
w_t_coeffs = q * (2 * ts + omega * sin(omega * ts));
w_tt_coeffs = q * (2 + omega^2 * cos(omega * ts));

% Define substrate coefficients such that w = a + b * x^2
as = w_coeffs * L^2;
a_ts = w_t_coeffs * L^2;
a_tts = w_tt_coeffs * L^2;

bs = - as / L^2;
b_ts = - a_ts / L^2;
b_tts = - a_tts / L^2;


%% Load in substrate coefficients
SubstrateCoefficients ...
    = substratecoefficients(as, bs, a_ts, b_ts, a_tts, b_tts, epsilon);
TimeDependents = timedependents(ts, SubstrateCoefficients);

%% Load in zero substrate coefficients
zeroTerm = zeros(size(ts));
StationarySubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm, zeroTerm, ...
    zeroTerm, zeroTerm, epsilon);
StationaryTimeDependents = timedependents(ts, StationarySubstrateCoefficients);

%% Load in forces
[Fs_composite, Fs_outer, Fs_inner] ...
    = substrateforce(ts, TimeDependents, epsilon);
[Stationary_Fs_composite, Stationary_Fs_outer, Stationary_Fs_inner] ...
    = substrateforce(ts, StationaryTimeDependents, epsilon);

%% Initialise figNo
figNo = 1;

%% Plot substrate
figure(figNo);
figNo = figNo + 1;
plot(ts, w_coeffs);
title("Substrate");

%% Plot turnover point
figure(figNo);
figNo = figNo + 1;
ds = TimeDependents.ds;
dStat = StationaryTimeDependents.ds;
hold on;
plot(ts, dStat);
plot(ts, ds);
legend(["Stationary", "Moving"]);
title("Turnover point");


%% Plot turnover derivative
figure(figNo);
figNo = figNo + 1;
d_ts = TimeDependents.d_ts;
d_tStat = StationaryTimeDependents.d_ts;
hold on;
plot(ts, d_tStat);
plot(ts, d_ts);
legend(["Stationary", "Moving"]);
title("Turnover point derivative");

%% Plot jet thickness point
figure(figNo);
figNo = figNo + 1;
Js = TimeDependents.Js;
JStat = StationaryTimeDependents.Js;
hold on;
plot(ts, JStat);
plot(ts, Js);
legend(["Stationary", "Moving"]);
title("Jet thickness");

%% Plot outer force
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, Stationary_Fs_outer);
plot(ts, Fs_outer);
legend(["Stationary", "Moving"]);
title("Outer force");

%% Plot composite force
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, Stationary_Fs_composite);
plot(ts, Fs_composite);
legend(["Stationary", "Moving"]);
title("Composite force");