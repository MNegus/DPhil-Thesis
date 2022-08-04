%% save_plate_turnovers.m
% Saves the turnover points for the validation cases with the plate code.
function save_plate_turnovers(master_dir)
    close all;

    %% Find a couple of examples
    max_levels = 10 : 14;
%     master_dir = "/scratch/negus/stationary_plate_maxlevel_validation";

    %% Options
    output_range = 0 : 800;
    Deltas = 6 ./ 2.^max_levels;
    dt = 1e-3;
    transpose_coordinates = false;
    shift = 0;
    % curvature_time = max(output_range) * dt;
    curvature_time = 0.45;

    %% Loop over both the 2D and the axisymmetric cases
    for axi = [0, 1]
        % Name of parent directory
        parent_dir = sprintf("%s/axi_%d", master_dir, axi);

        %% Loop over all the levels
        for m = 1 : length(max_levels)
            max_level = max_levels(m)
            Delta = Deltas(m);

            % Directories for this level
            level_dir = sprintf("%s/max_level_%d", parent_dir, max_level);
            interface_dir = append(level_dir, "/interfaces");

            % Find the turnover points
            [ts, ds, Js] = find_turnover_points(interface_dir, ...
                output_range, Delta, dt, transpose_coordinates, shift, ...
                curvature_time);

            % Save the turnover points as a text file
            combined = [ts; ds; Js];
            writematrix(combined', sprintf("%s/turnover_points.csv", level_dir));
        end
    end
end