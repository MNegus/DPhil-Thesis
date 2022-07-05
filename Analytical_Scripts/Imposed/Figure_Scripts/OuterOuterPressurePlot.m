function OuterOuterPressurePlot(dimension)

%% OuterPressure2D
% Script to produce a figure of the streamlines and pressure in the outer
% region for the 2D impact case. Appears in Section 3.3.6, under the Wagner
% theory chapter. 

close all;

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(groot, 'DefaultLegendInterpreter', 'latex');

%% Load in color map
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;
colormap(cmap);

%% Parameters for each dimension
if dimension == "2D"
    dirName = "Two-dimensional_Figures";
    minP = 0;
    maxP = 250;
else
    dirName = "Axisymmetric_Figures";
    minP = 0;
    maxP = 5000;
end

%% Parameter definitions
t = 0.5; % Times

%% Load in substrate functions
if dimension == "2D"
    substrateType = "curved";
    Pfun = @(A, X, Z) Pfun_2D(A, X, Z);
else
    substrateType = "flat";
    Pfun = @(A, X, Z) Pfun_axi(A, X, Z);
end
StationaryFunctions = substratefunctions("stationary", dimension);
SubstrateFunctions = substratefunctions(substrateType, dimension);

%% Save A coefficients
A_Stationary = StationaryFunctions.A(t);
A_Substrate = SubstrateFunctions.A(t);

%% Plotting parameters
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

if dimension == "2D"
    cb.Label.String = '$P_0(x, z, t)$';
    cb.XTick = [0, 50, 100, 150, 200, 250];
else
    cb.Label.String = '$P_0(r, z, t)$'; 
    cb.XTick = [0, 2500, 5000];
end
    

%% Whole figure settings
grid on;
if dimension == "2D"
    xlabel("$x$");
else
    xlabel("$r$");
end
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
filename = sprintf("OuterOuterWhole_%s", dimension);
savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);


%% Wedge plot
plotLim = 0.1;
rMin = 1 - 2 * plotLim;
ang = asin(plotLim / rMin);
numRs = 5e2;
numThetas = 5e2;

rs = sqrt(linspace((1 - 2 * plotLim).^2, 1, numRs));
thetas = linspace(-ang, ang, numThetas);

levels = linspace(minP, maxP, 30);

figure(2);
hold on;
colormap(cmap);

% Stationary substrate pressure (left)
thetaLeft = linspace(-ang, 0, numThetas / 2);
[R, Theta] = meshgrid(rs, thetaLeft);
X = R .* sin(Theta);
Z = 1 - R .* cos(Theta);
Pleft = Pfun(A_Stationary, X, Z);
Pleft(Pleft < 0) = 0;
contourf(X, Z, Pleft, levels, 'Edgecolor', 'None');

% Moving substrate pressure (right)
thetaRight = linspace(0, ang, numThetas / 2);
[R, Theta] = meshgrid(rs, thetaRight);
X = R .* sin(Theta);
Z = 1 - R .* cos(Theta);
Pright = Pfun(A_Substrate, X, Z);
Pright(Pright < 0) = 0;
contourf(X, Z, Pright, levels, 'Edgecolor', 'None');

% Add black line for interface
plot(sin(thetas), 1 - cos(thetas), 'Color', 'Black');

% Add white dividing line
xline(0, 'color', 'white');

xlim([-plotLim, plotLim]);
ylim([0, 2 * plotLim]);
pbaspect([1 1 1]);
grid on;
if dimension == "2D"
    xlabel("$x$");
else
    xlabel("$r$");
end
ylabel("$z$");

set(gcf,'position', [0, 0, 300, 300]);

filename = sprintf("OuterOuterZoomed_%s", dimension);
savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);


%% Function definitions
function Ps = Pfun_2D(A, X, Z)
%PFUN Outputs the pressure 
    Ps = A * (1 - X.^2 - (Z - 1).^2) ./ (2 * (X.^2 + Z.^2));
end

function Ps = Pfun_axi(A, X, Z)
    Ps = A * (1 - X.^2 - (Z - 1).^2) ./ (2 * (X.^2 + Z.^2)).^(3/2);
end

end