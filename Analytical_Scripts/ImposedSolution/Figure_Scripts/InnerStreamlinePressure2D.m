%% InnerStreamlinePressure
% Script to produce a figure of the streamlines and pressure in the inner
% region for the 2D impact case. Appears in Section 3.3.4, under the Wagner
% theory chapter. 
% 
% For the figure, we neglect any motion of the membrane, taking w == 0. 
%
% ADD DESCRIPTION AS TO WHAT WE ARE DOING IN TERMS OF PLOTTING STREAMLINES
% ETC

clear;
close all;

%% Figure options
set(0,'defaultTextInterpreter','latex'); %trying to set the default
set(0,'defaultAxesFontSize', 18);
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(groot, 'DefaultLegendInterpreter', 'latex');

xMax = 1;
zMax = 1;

%% Load in color map
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;

%% Parameter definitions
% Time dependent quantities
d = 1; % Let the turnover point be at 1
t = d^2 / 4; % Solves for t as a function of d
d_t = 1 / sqrt(t); % Time derivative of d
J = pi * d / (8 * d_t^2); % Jet thickness

% Anonymous functions
tzeta = @(eta) (J / pi) * (eta + 4 * 1i * sqrt(eta) - log(eta) + 1i * pi - 1);

% Difference between kappa values
dkappa = 0.25 * pi;

%% Figure properties
figure1 = figure();
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Set colour map to cool-warm
colormap(cmap);

% Axes labels
xlabel('$\tilde{x}$');
ylabel('$\tilde{z}$');

% Show axes grid
% grid on;

% Create the etas figure
% figure(2);
% hold on;

%% Variables (found from running before)
% maxEtaReal = -1e9;
% maxEtaImag = -1e9;
% minEtaReal = 1e9;
% minEtaImag = 1e9;

maxEtaReal = 56;
minEtaReal = -24;
maxEtaImag = 30;
minEtaImag = 1e-22;

%% Plot the pressure distribution
% Determine rMax
rMax = max(abs([minEtaReal, minEtaImag, maxEtaReal, maxEtaImag]));

% Discretise rs, clustered at r = 0
rs = exp(linspace(-100, 0, 1e3)) * rMax;

% Discretise thetas from 0 to pi
thetas = linspace(0, pi, 1e3);

% Create a polar meshgrid
[R, Thetas] = meshgrid(rs, thetas);

% Determine eta values
Etas = R .* exp(1i * Thetas);

% Determine meshgrid values for TZetas
TZetas = tzeta(Etas);
X = real(TZetas);
Z = imag(TZetas);

% Determine pressure values (divided by d_t^2)
P = (1 / 2) * (1 - abs((1 + 1i * sqrt(Etas)) ./ (1 - 1i * sqrt(Etas))).^2);

% Pressure levels
noFillConts = 30; % Number of filled contours
% levels = exp(linspace(-6, 1.8, noFillConts));
levels = 40;

% Plot the pressure
% figure(1);
contourf(X / J, Z / J, P, levels, 'Edgecolor', 'None');

%% Plot the free-surface
noPoints = 1e3;
xiMax = 100;

% Lower values
xisLower = linspace(-xiMax, 0, noPoints); 
xsLower = (J / pi) * (exp(xisLower) - xisLower - 1);
hsLower = (J / pi) * (pi + 4 * exp(xisLower / 2));

% Upper values
xisUpper = linspace(0, xiMax, noPoints / 2); 
xsUpper = (J / pi) * (xisUpper - log(1 + xisUpper));
hsUpper = (J / pi) * (pi + 4 * sqrt(xisUpper + 1));

% Combine and plot
% figure(1);
plot([xsLower, xsUpper] / J, [hsLower, hsUpper] / J, 'color', 'black', 'linewidth', 2);

%% Plot the jet-streamlines
% kappas = linspace(0, pi, 10); % kappa from 0 to pi
kappas = dkappa : dkappa : pi - dkappa;
rLower = exp(xisLower); % r values
rMin = 0.1; 

% Loop over the values of kappa
for kappa = kappas
    
    %% Lower side - into the jet
    % Solve for thetaSwitch, the point at which r == 1
    thetaSwitch = fsolve(@(theta) kappa - theta - rMin * sin(theta), 0.5 * pi);
    
    % Solve for the lower theta values
    theta0 = kappa * ones(size(rLower));
    thetaLower = fsolve(@(theta) theta - kappa + rLower .* sin(theta), theta0);
    
    % Determine the lower etas
    etaLower = rLower .* exp(1i * thetaLower);
    
    %% Upper side - into the outer region
    % Linear discretisation for theta
    thetaUpper = linspace(0, thetaSwitch, noPoints / 2);
    
    % Determine r from thetaUpper
    rUpper = (kappa - thetaUpper) ./ sin(thetaUpper);
    
    % Determine the lower etas
    etaUpper = rUpper .* exp(1i * thetaUpper);
    
    %% Create the continuous discretisation for etas and tzetas
    etas = [etaLower, etaUpper];
    tzetas = tzeta(etas);
    
    %% Find x and z values
    xs = real(tzetas);
    zs = imag(tzetas);
    
    %% Sort values in increasing z order
    [zs, sortIdxs] = sort(zs);
    xs = xs(sortIdxs);
    etas = etas(sortIdxs);
    
    %% Restrict to be within the axes limits
    restrictIdxs = abs(xs) < 1.5 * xMax & zs < 1.5 * zMax;
    
    %% Update eta limits
%     maxEtaReal = max(max(real(etas(restrictIdxs))), maxEtaReal);
%     minEtaReal = min(min(real(etas(restrictIdxs))), minEtaReal);
%     maxEtaImag = max(max(imag(etas(restrictIdxs))), maxEtaImag);
%     minEtaImag = min(min(imag(etas(restrictIdxs))), minEtaImag);
    
    %% Plot the streamline
%     figure(1);
    plot(xs(restrictIdxs) / J, zs(restrictIdxs) / J, 'color', 0.25 * [1 1 1], 'Linewidth', 2);
    
    %% Plot etas
%     figure(2);
%     plot(real(etas(restrictIdxs)), imag(etas(restrictIdxs)));
end

%% Plot the stagnation line
kappa = pi; % kappa is pi on the stagnation line
thetas = linspace(0, pi); % Linearly distributing theta is okay now
rs = (kappa - thetas) ./ sin(thetas); % Solution for r
etas = rs .* exp(1i * thetas); % Solution for eta

% Solution for tzeta, where we manually set the final point to be at the
% stagnation point
tzetas = tzeta(etas);
tzetas(end) = tzeta(-1);

% Plot the stagnation line
% figure(1);
plot(real(tzetas) / J, imag(tzetas) / J, 'color', 'black', 'linestyle', '--', 'linewidth', 2);

%% Plot the outer-streamlines
% MAKE THE DIFFERENCE BETWEEN EACH KAPPA THE SAME
kappas = pi + 0.1 * dkappa : dkappa : 10 * pi
thetas = linspace(0, pi); % Linearly distributing theta is okay now

for kappa = kappas(2 : end)
    rs = (kappa - thetas) ./ sin(thetas); % Solution for r
    etas = rs .* exp(1i * thetas); % Solution for eta

    % Solution for tzeta
    tzetas = tzeta(etas);

    % Find x and z values
    xs = real(tzetas);
    zs = imag(tzetas);
    
    % Restrict to be within the axes limits
    restrictIdxs = abs(xs) < 1.5 * xMax & zs < 1.5 * zMax;
    
%     % Update eta limits
%     maxEtaReal = max(max(real(etas(restrictIdxs))), maxEtaReal);
%     minEtaReal = min(min(real(etas(restrictIdxs))), minEtaReal);
%     maxEtaImag = max(max(imag(etas(restrictIdxs))), maxEtaImag);
%     minEtaImag = min(min(imag(etas(restrictIdxs))), minEtaImag);
    
    % Plot the streamline
%     figure(1);
    plot(xs(restrictIdxs) / J, zs(restrictIdxs) / J, 'color', 0.25 * [1 1 1], 'Linewidth', 2);
    
    % Plot etas
%     figure(2);
%     plot(real(etas(restrictIdxs)), imag(etas(restrictIdxs)));
    
end

%% Plot colour bar for the pressure plot
cb = colorbar('Location', 'Northoutside');
cb.Label.String = '$\tilde{p}_0(\tilde{x}, \tilde{z}, t) / \dot{d}_0(t)^2$';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'latex';

% Set position of colour bar to be above the axes at the centre (a bit
% hacky but works for now)
x1 = get(gca,'position');
x = get(cb,'Position');
x(3) = 0.5 * x(3);
x(2) = 0.99 * x(2);
x(1) = 2.475 * x(1);
set(cb,'Position',x)
set(gca,'position',x1)



%% Figure scalings
figure(1);
xlim([-xMax / J, xMax / J]);
ylim([0, zMax / J]);
pbaspect([1 zMax / (2 * xMax) 1]);

% Set tick values
xTicks = ["$-8J(t)$", "$-4J(t)$", "$0$", "$4J(t)$", "$8J(t)$"];
yTicks = ["$0$", "$2J(t)$", "$4J(t)$", "$6J(t)$", "$8J(t)$", "$10J(t)$"];
set(gca,'xtick', [-8 : 4 : 8],'xticklabel', xTicks);
set(gca,'ytick', [0 : 2 : 10],'yticklabel', yTicks);
% Add J(t) labellings to x and y limits


% Set figure size
set(gcf,'position', [0, 0, 1200, 600]);
box on;

% Set rendered to Painters (incredibly slow but makes the figures better)
% set(gcf, 'Renderer', 'Painters');


%% Create figures
filename = "InnerStreamlinePressure2D";
dirName = "Two-dimensional_Figures";
savefig(gcf, sprintf("%s/fig/%s.fig", dirName, filename));
exportgraphics(gcf, sprintf("%s/png/%s.png", dirName, filename), 'Resolution', 300);
exportgraphics(gcf,sprintf("%s/eps/%s.eps", dirName, filename), 'Resolution', 300);


