%% mass_validation.m
% Validates the mass of the droplet in the imposed deformable substrate
% case

clear;

% Load in red-blue colour map
cmap_mat = matfile('red_blue_cmap.mat');
cmap = cmap_mat.cmap;

fontsize = 22;

%% Parameters
L = 2;
DELTA_T = 1e-4;
IMPACT_TIME = 0.125;
T_MAX = 0.8;
ts = - IMPACT_TIME : DELTA_T : T_MAX - IMPACT_TIME;
ts_analytical = 0 : DELTA_T : T_MAX - IMPACT_TIME;
MAX_TIMESTEP = T_MAX / DELTA_T;
NO_TIMESTEPS = length(ts);
freq = 100; % Frequency

%% Directory names
master_dir = "/scratch/negus/DNS_Chapter_validation";
type = ["imposed_membrane_0.0025"];
levels = 10 : 14;

%% Set up colors
color_idxs = floor(linspace(1, length(cmap), length(levels)));
no_levels = length(levels);
colors = ones(no_levels - 1, 3);

for q = 1 : no_levels - 1
    colors(q, :) = cmap(color_idxs(q), :);
end

%% Analytical solution
epsilon = 1; % Analytical parameter for small time 
a = 0.0025; % Imposed parameter
[ws, w_ts, w_tts] = ImposedPlate(ts_analytical, a);
lambda = pi / (2 * L);

% Solution for ds
options = optimoptions('fsolve', 'TolFun', 1e-10, 'TolX', 1e-10);
dsGuess = 2 * sqrt(ts_analytical);
ds_zero_fun = @(d) ts_analytical - d.^2 / 4 - ws .* besselj(0, epsilon * lambda * d);
ds = fsolve(ds_zero_fun, dsGuess, options);

% Solution for fluxes
fluxes = 2 * w_ts .* sin(lambda * ds) / lambda;

% Solution for mass loss
massLoss = -cumtrapz(ts_analytical, fluxes);
AsAnalytical = massLoss / pi;

%% Plot all cases
close all;

imposed_coeffs = ["0"; "0.0025"];
for imposedIdx = 1 : 2
    
    imposed_coeff = imposed_coeffs(imposedIdx);
    
    parent_dir = sprintf("%s/%s_maxlevel_validation/imposed_coeff_%s", ...
        master_dir, type, imposed_coeff);

    % Create figure
    close(figure(1));
    figure(1);
%     figure(imposedIdx);
    hold on;

    % Matrix to hold the turnover points
    dns_masses = zeros(NO_TIMESTEPS, length(levels));

    % Loop over the levels
    for level_idx = 1 : length(levels)
        level = levels(level_idx);
        level_dir = sprintf("%s/max_level_%d", parent_dir, level);

        mass_mat = readmatrix(sprintf("%s/output.txt", level_dir));
        ts = mass_mat(1 : NO_TIMESTEPS, 1) - IMPACT_TIME;
        mass = 2 * mass_mat(1 : NO_TIMESTEPS, 2) / (2 * pi);
        
        mass_ratio = (mass - pi) / pi;

        dns_masses(:, level_idx) = mass_ratio;

        scatter(ts(1 : freq : end), mass_ratio(1 : freq : end), [], cmap(color_idxs(level_idx), :), ...
            'linewidth', 2, 'Displayname', sprintf("$m$ = %d", level));
        drawnow;
    end
    % Analytical solution
    if imposed_coeff == '0'
        plot(ts_analytical, zeros(size(ts_analytical)), ...
            'linestyle', '--', 'linewidth', 2, 'color', 'black', ...
            'displayname', 'Analytical');
    else
        plot(ts_analytical, AsAnalytical, 'linestyle', '--', 'linewidth', 2, ...
            'color', 'black', 'displayname', 'Analytical');
    end

    grid on;
    xlim([-0.32, 0.8]);
    ylim([-4, 5] * 10^(-3));
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    xlabel("$t$", 'interpreter', 'latex', "Fontsize", fontsize);
    ylabel("$R_m(t)$", 'interpreter', 'latex');
    legend('location', 'southwest', 'interpreter', 'latex');

    %% Plot L2 norm error
    axes('Position',[.36 .6 .3 .3])
    box on
    hold on;
    mass_max = dns_masses(:, length(levels));


    norms = zeros(length(levels) - 1, 1);
    for level_idx = 1 : length(levels) - 1
        level = levels(level_idx);

        % Determines L2-norm difference
        diff = dns_masses(:, level_idx) - mass_max;
        norms(level_idx) = sqrt(sum(diff.^2) / length(diff));
    end
    plot(levels(1 : end - 1), norms, 'color', 'black', 'linewidth', 2);

    sz = 100;
    scatter(levels(1 : end - 1), norms, sz, colors, 'filled');

    set(gca, 'yscale', 'log');
    xlim([min(levels) - 0.5, max(levels) - 0.5])
    ylim([10^-5.5, 10^-2.5]);
    xticks(levels(1 : end - 1));
    yticks([10^-5, 10^-4, 10^-3]);

    grid on;
    set(gca, "ticklabelinterpreter", "latex", "Fontsize", fontsize);
    xlabel("$m$", 'interpreter', 'latex');
    ylabel("$||R_m - R_{14}||_2$", 'interpreter', 'latex');

    %% Overal figure settings
    x0=400;
    y0=400;
    height=800;
    width=600;

    set(gcf,'position',[x0,y0,width,height]);
    

    if imposed_coeff == "0"
        figname = "dns_validation_figures/areas/DNSAreas_Stationary";
    else
        figname = "dns_validation_figures/areas/DNSAreas_Imposed";
    end
    set(gcf, 'Renderer', 'Painters');
    pause(1.5);
    exportgraphics(gcf, sprintf("%s.png", figname), "Resolution", 300);
    savefig(gcf, sprintf("%s.fig", figname));

end

    