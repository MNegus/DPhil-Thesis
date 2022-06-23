clear;
close all;

% addpath("Pressures");
addpath("Forces");

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');


%% Parameters
tmax = 1;
ts = linspace(0, tmax, 1e3)';

%% Substrate parameters, 
[epsilon, L, q, omega] = plateparameters(); % Substrate parameters

[ws, w_ts, w_tts] = flatsubstrate(ts, q, omega);


%% Load in substrate coefficients
SubstrateCoefficients = substratecoefficients(ws, w_ts, w_tts);
TimeDependents = timedependents(ts, SubstrateCoefficients);

%% Load in zero substrate coefficients
zeroTerm = zeros(size(ts));
StationarySubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm);
StationaryTimeDependents = timedependents(ts, StationarySubstrateCoefficients);

%% Load in forces
[Fs_composite, Fs_outer, Fs_inner] ...
    = substrateforce(ts, TimeDependents, epsilon);
[Stationary_Fs_composite, Stationary_Fs_outer, Stationary_Fs_inner] ...
    = substrateforce(ts, StationaryTimeDependents, epsilon);

%% Initialise figNo
figNo = 1;



%% Plot forces
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, Stationary_Fs_composite);
plot(ts, Fs_composite);

