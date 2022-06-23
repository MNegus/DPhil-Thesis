clear;
close all;

%% Load in coarse color map
mapObj = load("red_blue_cmap.mat");
coarsecmap = mapObj.cmap;
coarseIdxs = 1 : length(coarsecmap);

%% Smoothen color map using spline interpolation
noPoints = 1e3;
fineIdxs = linspace(1, length(coarsecmap), noPoints);
cmap = zeros(noPoints, 3);

for col = 1 : 3
    % Smoothen column using spline
    finePoints = interp1(coarseIdxs, coarsecmap(:, col), fineIdxs, 'spline');
    
    % Ensure all points are between 0 and 1
    finePoints(finePoints > 1) = 1;
    finePoints(finePoints < 0) = 0;
    
    % Save column
    cmap(:, col) = finePoints;
end

% Save matrix file
save('fine_red_blue_cmap.mat', 'cmap');
