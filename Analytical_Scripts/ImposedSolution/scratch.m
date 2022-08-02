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
tmax = 1;
ts = linspace(0, tmax, 1e3)';
[epsilon, k, q, omega] = substrateparameters();

%% Load in substrate functions
StationaryFunctions = substratefunctions("stationary");
FlatFunctions = substratefunctions("flat");
CurvedFunctions = substratefunctions("curved");

%% Substrate parameters, such that w(x, t) = w(t) * (L^2 - x^2);
% [epsilon, L, q, omega] = quadraticparameters(); % Substrate parameters
% 
% [aQuads, a_tQuads, a_ttQuads, bQuads, b_tQuads, b_ttQuads] ...
%     = quadraticsubstrate(ts, L, q, omega); % Quadratic substrate coefficients
% 
% [aFlats, a_tFlats, a_ttFlats, bFlats, b_tFlats, b_ttFlats] ...
%     = flatsubstrate(ts, q, omega); % Flat substrate coefficients


%% Load in substrate coefficients
% QuadSubstrateCoefficients ...
%     = substratecoefficients(aQuads, bQuads, a_tQuads, b_tQuads, a_ttQuads, b_ttQuads, epsilon);
% QuadTimeDependents = timedependents(ts, QuadSubstrateCoefficients);
% 
% FlatSubstrateCoefficients ...
%     = substratecoefficients(aFlats, bFlats, a_tFlats, b_tFlats, a_ttFlats, b_ttFlats, epsilon);
% FlatTimeDependents = timedependents(ts, FlatSubstrateCoefficients);

%% Load in zero substrate coefficients
% zeroTerm = zeros(size(ts));
% StationarySubstrateCoefficients ...
%     = substratecoefficients(zeroTerm, zeroTerm, zeroTerm, zeroTerm, ...
%     zeroTerm, zeroTerm, epsilon);
% StationaryTimeDependents = timedependents(ts, StationarySubstrateCoefficients);

%% Load in forces
[Stationary_Fs_composite, Stationary_Fs_outer, Stationary_Fs_inner] ...
    = substrateforce(ts, StationaryFunctions, epsilon);
[Flat_Fs_composite, Flat_Fs_outer, Flat_Fs_inner] ...
    = substrateforce(ts, FlatFunctions, epsilon);
[Quad_Fs_composite, Quad_Fs_outer, Quad_Fs_inner] ...
    = substrateforce(ts, CurvedFunctions, epsilon);


%% Initialise figNo
figNo = 1;

%% Plot substrate
% figure(figNo);
% figNo = figNo + 1;
% plot(ts, aQuads);
% title("Substrate");

%% Plot substrate velocity
% figure(figNo);
% figNo = figNo + 1;
% plot(ts, a_tQuads);
% title("Substrate velocity");

%% Plot substrate velocity
% figure(figNo);
% figNo = figNo + 1;
% plot(ts, a_ttQuads);
% title("Substrate acceleration");

%% Plot turnover point
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryFunctions.d(ts));
plot(ts, FlatFunctions.d(ts));
plot(ts, CurvedFunctions.d(ts));
legend(["Stationary", "Flat", "Quadratic"]);
title("Turnover point");


%% Plot turnover derivative
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryFunctions.d_t(ts));
plot(ts, FlatFunctions.d_t(ts));
plot(ts, CurvedFunctions.d_t(ts));
legend(["Stationary", "Flat", "Quadratic"]);
title("Turnover point derivative");

%% Plot jet thickness point
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryFunctions.J(ts));
plot(ts, FlatFunctions.J(ts));
plot(ts, CurvedFunctions.J(ts));
legend(["Stationary", "Flat", "Quadratic"]);
title("Jet thickness");

%% Plot outer force
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, Stationary_Fs_outer);
plot(ts, Flat_Fs_outer);
plot(ts, Quad_Fs_outer);
legend(["Stationary", "Flat", "Quadratic"]);
title("Outer force");

%% Plot composite force
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, Stationary_Fs_composite);
plot(ts, Flat_Fs_composite);
plot(ts, Quad_Fs_composite);
legend(["Stationary", "Flat", "Quadratic"]);
title("Composite force");

%% Plot A values
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryFunctions.A(ts));
plot(ts, FlatFunctions.A(ts));
plot(ts, CurvedFunctions.A(ts));
legend(["Stationary", "Flat", "Quadratic"]);
title("A(t)");

%% Plot B values
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryFunctions.B(ts));
plot(ts, FlatFunctions.B(ts));
plot(ts, CurvedFunctions.B(ts));
legend(["Stationary", "Flat", "Quadratic"]);
title("B(t)");

%% Plot C values
figure(figNo);
figNo = figNo + 1;
hold on;
plot(ts, StationaryFunctions.C(ts));
plot(ts, FlatFunctions.C(ts));
plot(ts, CurvedFunctions.C(ts));
legend(["Stationary", "Flat", "Quadratic"]);
title("C(t)");
