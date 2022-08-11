%% save_membrane_turnovers.m
% Saves the turnover points for the validation cases with the membrane code.

function save_membrane_turnovers(master_dir)
    close all;

    %% Find a couple of examples
    max_levels = 10 : 14;

    %% Options
    output_range = 0 : 10 : 8000;
    Deltas = 6 ./ 2.^max_levels;
    dt = 1e-4;
    transpose_coordinates = false;
    shift = 0;
    curvature_time = 0.45;


    %% Loop over both imposed cases
    % for acc = ["0.0"; "0.05"]
    accs = ["0", "0.0025"];
    for k = 1 : 2
        % Name of parent directory
        acc = accs(k);
        parent_dir = sprintf("%s/imposed_coeff_%s", master_dir, acc);

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