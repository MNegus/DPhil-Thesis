clear;
close all;

addpath("Pressures");
addpath("Forces");

%% Parameters
epsilon = 0.1;
% tmax = 1 / (4 * epsilon^2);
tmax = 1;
ts = linspace(1e-10, tmax, 1e3)';

%% Load in plate positions and time-dependents
ws = zeros(size(ts));
w_ts = zeros(size(ts));
w_tts = zeros(size(ts));

PlatePositions = platepositions(ws, w_ts, w_tts);
TimeDependents = timedependents(ts, PlatePositions);

%% Determine forces
[Fs_composite, Fs_outer, Fs_inner] ...
    = substrateforce(PlatePositions, TimeDependents, epsilon);

%% Plot forces
figure();
hold on;
plot(ts, Fs_outer);
plot(ts, Fs_composite);


