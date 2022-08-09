%% OuterOuterPressureComparison.m
%

clear;
close all;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/");
addpath("../Analytical_Scripts/PlateSolution/");
addpath("../Analytical_Scripts/Pressures/");
addpath("../Analytical_Scripts/Forces/");

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;
colormap(cmap);

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

%% Parameters
% Plate parameters
ALPHA = 2;
BETA = 0;
GAMMA = 20;
L = 2; 

% Computational parameters
DELTA_T = 1e-3;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
MAX_TIMESTEP = T_MAX / DELTA_T;

% DNS directories
dns_master_dir = "/home/michael/scratch/DPhil_DNS_Data/";
stat_dir = sprintf("%s/Stationary_Plate/axi", dns_master_dir);
moving_dir = sprintf("%s/Moving_Plate/ALPHA-%g_BETA-%g_GAMMA-%g", dns_master_dir, ...
    ALPHA, BETA, GAMMA);

% Analytical parameters
epsilon = 1;

% Analytical xs array
xsAnalytical = linspace(0, L, 1e3);

%% Load stationary analytical solutions
tMax = 1 / 3; % Max such that 1 = d(t)
tsStat = linspace(0, tMax, 1e3);

% Find substrate functions
StatSubstrateFunctions = platesubstratefunctions(tsStat, ...
    zeros(size(tsStat)), zeros(size(tsStat)), zeros(size(tsStat)), epsilon);

% Find ds
dsStat = StatSubstrateFunctions.d(tsStat);

%% Load moving analytical solutions
tMax = 1.5 * (1 / 3); % A bit after the stationary max value

% Solve plate equation
[tsComp, wsComp, w_tsComp, w_ttsComp] ...
    = PlateSolution(tMax, ALPHA, BETA, GAMMA, epsilon, "composite");

% Find substrate functions
MovingSubstrateFunctions = platesubstratefunctions(tsComp, ...
    wsComp, w_tsComp, w_ttsComp, epsilon);

% Find where turnover point reaches 1
dsComp = MovingSubstrateFunctions.d(tsComp);
tIdxMaxComp = sum(dsComp <= 1);

% Restrict solutions temporally
tsComp = tsComp(1 : tIdxMaxComp);
wsComp = wsComp(1 : tIdxMaxComp);
w_tsComp = w_tsComp(1 : tIdxMaxComp);
w_ttsComp = w_ttsComp(1 : tIdxMaxComp);
dsComp = dsComp(1 : tIdxMaxComp);

%% Select a time
t = 0.075;

%% Save A coefficients
A_Stationary = StatSubstrateFunctions.A(t);
A_Substrate = MovingSubstrateFunctions.A(t);

%% Plotting parameters
minP = 0;
maxP = 4.8823;
levels = linspace(minP, maxP, 100);

%% Surf plot for pressure
numRs = 1e2;
numThetas = 1e3;
rs = sqrt(sqrt(linspace(0, 1, numRs)));
thetas = linspace(0, 2 * pi, numThetas);

figure(1);
hold on;

% Stationary substrate pressure (left)
thetaLeft = linspace(pi, 2 * pi, numThetas / 2);
[R, Theta] = meshgrid(rs, thetaLeft);
X = R .* sin(Theta);
Z = 1 - R .* cos(Theta);
Pleft = Pfun(A_Stationary, X, Z);
Pleft(Pleft < 0) = 0;
contourf(X, Z, Pleft, levels, 'Edgecolor', 'None');

% Moving substrate pressure (right)
thetaRight = linspace(0, pi, numThetas / 2);
[R, Theta] = meshgrid(rs, thetaRight);
X = R .* sin(Theta);
Z = 1 - R .* cos(Theta);
Pright = Pfun(A_Substrate, X, Z);
Pright(Pright < 0) = 0;
contourf(X, Z, Pright, levels, 'Edgecolor', 'None');

plot(sin(thetas), 1 + cos(thetas), 'Color', 'Black');

xlim([-1, 1]);
ylim([0, 2]);
pbaspect([1 1 1]);

% Colourbar 
cb = colorbar('Location', 'Northoutside');
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';
% 
caxis([minP maxP]);

% if dimension == "2D"
%     cb.Label.String = '$P_0(x, z, t)$';
% %     cb.XTick = [0, 50, 100, 150, 200, 250];
% else
%     cb.Label.String = '$\epsilon P_0(r, z, t)$'; 
% %     cb.XTick = [0, 2500, 5000] * epsilon;
% end
%     cb.XTick = [0, 50, 100, 150, 200, 250];

%% Whole figure settings
grid on;
xlabel("$r$");
ylabel("$z$");
set(gcf,'position', [0, 0, 800, 600]);

% Set position of colour bar to be above the axes at the centre (a bit
% hacky but works for now)
axesPos = get(gca,'position');
figPos = get(gcf, 'Position');
cbPos = get(cb,'Position');

cbPos(3) = 0.75 * cbPos(3); % Half the width of the colour bar
cbPos(1) = (1.3 * axesPos(3) / 8) + 0.5 * cbPos(3); % Set left of the colour bar to be 1/8 along the axis
cbPos(2) = 0.83;

set(cb,'Position',cbPos)
set(gca,'Position',axesPos);

% Axes limits
pertb = 0.1;
extent = 1 + pertb;

% Add white dividing line
xline(0, 'color', 'white');

% Axes limits
xlim([-extent, extent]);
ylim([0, 2 * extent]);

% Set position of colour bar to be above the axes at the centre (a bit
% hacky but works for now)
axesPos = get(gca,'position');
figPos = get(gcf, 'Position');
cbPos = get(cb,'Position');

cbPos(3) = 0.75 * cbPos(3); % Half the width of the colour bar
cbPos(1) = (2.25 * axesPos(3) / 8) + 0.5 * cbPos(3); % Set left of the colour bar to be 1/8 along the axis
cbPos(2) = 0.83;

set(cb,'Position',cbPos)
set(gca,'Position',axesPos);

% Export figure
% filename = sprintf("OuterOuterWhole_%s", dimension);
% savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
% exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
% exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);


%% Pressure function
function Ps = Pfun(A, X, Z)
    Ps = A * (1 - X.^2 - (Z - 1).^2) ./ (2 * (X.^2 + Z.^2)).^(3/2);
end
