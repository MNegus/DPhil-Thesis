clear;
close all;

%% Load in color map
mapObj = load("fine_red_blue_cmap.mat");
cmap = mapObj.cmap;

%% Plot the red-blue bar
figure(1);
colormap(cmap);

cb = colorbar('Location', 'Northoutside');
set(gca,'Visible',false);
cb.Ticks = [];
exportgraphics(gcf, "RedBlueColourBar.png", "Resolution", 300);
% savefig(cb, 'colourbar.png')

%% Plot the white-blue bar
% Find midpoint of the red-blue bar
midIdx = floor(length(cmap) / 2);

% Cutoff half the map
cutoffcmap = cmap(1 : midIdx, :);

% Flip the order of the map
cutoffcmap = flip(cutoffcmap);

figure(2);
colormap(cutoffcmap);

cb = colorbar('Location', 'Northoutside');
set(gca,'Visible',false);
cb.Ticks = [];
exportgraphics(gcf, "WhiteBlueColourBar.png", "Resolution", 300);