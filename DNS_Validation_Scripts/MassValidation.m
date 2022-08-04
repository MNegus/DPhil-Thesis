%% mass_validation.m
% Validates the mass of the droplet in the imposed deformable substrate
% case

clear;
close all;

% Load in red-blue colour map
cmap_mat = matfile("../fine_red_blue_cmap.mat");
cmap = cmap_mat.cmap;

% Adds analytical scripts to path
addpath("../Analytical_Scripts/Imposed/");

%% Figure options
fontsize = 11;
lineWidth = 1.25;
set(0,'defaultTextInterpreter','latex');
set(0,'defaultAxesFontSize', fontsize);
set(0, 'DefaultTextFontSize', fontsize);
set(0,'defaultLegendFontSize', fontsize, 'DefaultLegendFontSizeMode','manual');
set(0, 'defaultAxesTickLabelInterpreter', 'latex');
set(0, 'defaultFigureRenderer', 'painters');
set(0, 'DefaultLegendInterpreter', 'latex');

%% Parameters
L = 2;
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);
freq = 100; % Frequency

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
type = ["imposed_0.05_omega_4"];
levels = 10 : 14;
imposed_coeffs = ["0"; "0.05"];

%% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels - 1, 3);

for q = 1 : no_levels - 1
    colors(q, :) = cmap(color_idxs(q), :);
end

%% Find analytical solution
[epsilon, q, omega, ~, L] = substrateparameters("curvedDNS");
SubstrateFunctions = substratefunctions("curvedDNS", "2D");

% Find analytical times
tsAnalyticalUnRestricted = (1e-9 : DELTA_T : T_MAX - IMPACT_TIME) / epsilon^2;

% Load analytical turnover points
dsAnalyticalUnRestricted = SubstrateFunctions.d(tsAnalyticalUnRestricted);

% Find times restricted to where where turnover point less that 1
% tsAnalytical = tsAnalyticalUnRestricted(epsilon * dsAnalyticalUnRestricted <= 1);
% ds = dsAnalyticalUnRestricted(epsilon * dsAnalyticalUnRestricted <= 1);
tsAnalytical = tsAnalyticalUnRestricted;
ds = dsAnalyticalUnRestricted;


% Find fluxes leaving boundary
fluxes = 2 * q * (2 * tsAnalytical + omega * sin(omega * tsAnalytical)) ... 
    .* ds .* (1 - ds.^2 / (3 * L^2));

% FLUX ASSUMING ALL PLATE IS WETTED
fluxesMax = 2 * q * (2 * tsAnalytical + omega * sin(omega * tsAnalytical)) ... 
    .* 2 * L / 3;


%% Plot all cases
% Set up tiled figure
tileFig = tiledlayout(1, 2, 'Tilespacing', 'compact', 'padding', 'compact');

% Norms matrix
norms = zeros(length(levels) - 1, 2);

for imposedIdx = 1 : 2
    
    imposed_coeff = imposed_coeffs(imposedIdx);
    
    parent_dir = sprintf("%s/%s_maxlevel_validation/imposed_coeff_%s", ...
        master_dir, type, imposed_coeff);

    tileFig; % Change to tiled figure

    % Save parent directory
    parent_dir = sprintf("%s/%s_maxlevel_validation/imposed_coeff_%s", ...
        master_dir, type, imposed_coeff);
    
    % Figure titles
    if imposed_coeff == "0"
        titleStr = "(a) Stationary substrate setup.";
    else
        titleStr = "(b) Linearised boundary setup.";
    end

    %% Plot DNS solutions
    % Matrix to hold the turnover points
    dns_masses = zeros(NO_TIMESTEPS, length(levels));

    nexttile(imposedIdx);
    hold on;
    
    % Loop over the levels
    for level_idx = 1 : length(levels)
        level = levels(level_idx);
        level_dir = sprintf("%s/max_level_%d", parent_dir, level);

        mass_mat = readmatrix(sprintf("%s/output.txt", level_dir));
        ts = mass_mat(1 : NO_TIMESTEPS, 1) - IMPACT_TIME;
        mass = 2 * mass_mat(1 : NO_TIMESTEPS, 2) / (2 * pi);
        
        mass_ratio = (mass - pi) / pi;

        dns_masses(:, level_idx) = mass_ratio;
    
        sz = 15;
        h(level_idx) = scatter(ts(1 : freq : end), mass_ratio(1 : freq : end), sz, cmap(color_idxs(level_idx), :), ...
            'linewidth', lineWidth, 'Displayname', sprintf("$m$ = %d", level));
    end
    
    
    %% Determine L2 norms
    mass_max = dns_masses(:, length(levels)); 
    for level_idx = 1 : length(levels) - 1
        level = levels(level_idx);

        % Determines L2-norm difference
        diff = dns_masses(:, level_idx) - mass_max;
        norms(level_idx, imposedIdx) = sqrt(sum(diff.^2) / length(diff));
    end
    
    
    %% Plot analytical solution
    if imposed_coeff == "0"
        RsAnalytical = zeros(size(tsAnalytical));
        RsMax = zeros(size(tsAnalytical));
    else
        RsAnalytical = -cumtrapz(tsAnalytical, fluxes) / pi;
        RsMax = -cumtrapz(tsAnalytical, fluxesMax) / pi;
    end
    
    h(length(levels) + 1) = scatter(10, 10, [], 'white', 'Displayname', '');
        
    h(length(levels) + 2) = plot(tsAnalytical, RsAnalytical, 'linestyle', '--', 'linewidth', 2, ...
        'color', 'black', 'displayname', 'Analytical');

    h(length(levels) + 3) = plot(tsAnalytical, RsMax, 'linestyle', ':', 'linewidth', 2, ...
        'color', 'black', 'displayname', 'Maximum loss');

    %% Tile figure options
    grid on;
    xlim([-0.32, 0.8]);
    if imposedIdx == 1
        ylim([-2.5*10^-3, 0.5 * 10^-3]);
    else
        ylim([-0.1, 0.01]);
        yticks(-0.1 : 0.02 : 0); 
    end
%     ylim([-4, 5] * 10^(-3));
    
    % Axes labels
    xlabel("$t$");
    ylabel("$R_m(t)$");

    % Create title
    title({titleStr ''}, 'Fontsize', fontsize);
    set(gca, 'TitleFontSizeMultiplier', 1);
end

%% Set options for entire figure
% Set legend
lh = legend(h(1 : length(levels) + 3), 'Numcolumns', 3);
lh.Layout.Tile = 'South'; 

% Set size in inches
width = 6;
height = 4;
set(gcf,'units', 'inches', ...
    'position',[0.5 * width, 0.5 * height, width, height]);

%% Plot L2 norm error

% Loop over types
for imposedIdx = 1 : 2
    
    % Define axes options for the inset
    insetWidth = 0.12;
    insetHeight = 0.14;
    verticalPos = 0.4;
    if imposedIdx == 1
        horizontalPos = 0.19;
    else
        horizontalPos = 0.67;
    end

    % Define inset axes
    axes('Position',[horizontalPos, verticalPos, insetWidth, insetHeight]);
    box on;
    grid on;
    hold on;

    plot(levels(1 : end - 1), norms(:, imposedIdx), 'color', 'black', ...
        'linewidth', lineWidth);

    sz = 35;
    scatter(levels(1 : end - 1), norms(:, imposedIdx), sz, colors, 'filled');

    set(gca, 'yscale', 'log');
    xlim([min(levels) - 0.5, max(levels) - 0.5])
    ylim([10^-5.5, 10^-2.5]);
    xticks(levels(1 : end - 1));
    yticks([10^-5, 10^-4, 10^-3]);

    grid on;
    xlabel("$m$");
    ylabel("$||R_m - R_{14}||_2$");
end

%% Exporting figure

set(gcf, 'Renderer', 'Painters');
pause(1.5);

figname = "dns_validation_figures/areas/DNSAreas";
exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
savefig(gcf, sprintf("%s.fig", figname));


    