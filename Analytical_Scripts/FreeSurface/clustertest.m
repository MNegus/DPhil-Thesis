
%% Cluster at one end point
xMin = 0;
xMax = 1;
lambda = 10;
sigmas = linspace(1, 0, 1e3);
b = (xMax - xMin) / (exp(-lambda) - 1);
a = xMin - b;

xs = a + b * exp(-lambda * sigmas)

plot(xs, zeros(size(xs)), '-o')

%% Cluster at both end points
% xMin = 0.1;
% xMax = 1;
% lambda = 10;
% sigmas = linspace(-1, 1, 1e3);
% a = 0.5 * (xMin + xMax);
% b = 0.5 * coth(lambda) * (xMax - xMin);
% 
% xs = a + b * tanh(lambda * sigmas);
% plot(xs, zeros(size(xs)), '-o');
