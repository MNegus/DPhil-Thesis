%% OuterPressurePlot

clear;
close all;


addpath("OuterSolution");

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Load in color map
% mapObj = load("red_blue_cmap.mat");
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;

%% Parameters
tmax = 1;
ts = tmax;

%% Substrate parameters, 
[epsilon, L, q, omega] = plateparameters(); % Substrate parameters

[ws, w_ts, w_tts] = flatsubstrate(ts, q, omega);


%% Load in zero substrate coefficients
zeroTerm = zeros(size(ts));
SubstrateCoefficients ...
    = substratecoefficients(zeroTerm, zeroTerm, zeroTerm);
TimeDependents = timedependents(ts, SubstrateCoefficients);

%% Load in d, d_t, d_tt
d = TimeDependents.ds
d_t = TimeDependents.d_ts
d_tt = TimeDependents.d_tts

%% Plotting parameters

% Spatial quantities
noPoints = 1e3; % Number of grid points per dimension
rMax = 2 * d; % Maximum value of x
zMax = rMax; % Maximum value of z

rs = linspace(0, rMax, noPoints); % Array of points in x
zs = linspace(0, zMax, noPoints); % Array of points in z

% Filled contour options
noFillConts = 30; % Number of filled contours
scaleFun = @(P) log(P); % Scale function for the filled contours
inverseScaleFun = @(val) exp(val); % Inverse of scale function

%% Figure properties
% figure1 = figure();
% axes1 = axes('Parent',figure1);
% hold(axes1,'on');
figure(2);
hold on;

% Set colour map to cool-warm
colormap(cmap);

% Axes labels
xlabel('$\hat{r}$');
ylabel('$\hat{z}$');

% Show axes grid
% grid on;

%% Creates mesh grid of values
[R, Z] = meshgrid(rs, zs);
P = zeros(size(R));

sizeR = size(R);
for j = 1 : sizeR(1)
    for k = 1 : sizeR(2)
        r = R(j, k);
        z = Z(j, k);
        
        P(j, k) = outerpressure(r, z, d, d_t, d_tt);
        [j, k]
    end
end

%% 
close all;
figure(1);
colormap(cmap);
contourf(R, Z, P, 50, 'edgecolor', 'none');
