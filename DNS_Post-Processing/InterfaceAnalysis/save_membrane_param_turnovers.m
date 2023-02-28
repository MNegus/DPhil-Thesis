function save_membrane_param_turnovers(master_dir)
    close all;

    %% Options
    max_level = 14;
    output_range = 0 : 700;
    box_width = 24;
    Delta = box_width ./ 2^max_level;
    dt = 1e-3;
    transpose_coordinates = false;
    shift = 0;
    curvature_time = 0.45;

    %% Loop over subdirectories
    dirStruct = dir(master_dir);

    for k = 3 : length(dirStruct)
        param_dir = append(master_dir, "/", dirStruct(k).name)

        % Interface directory
        interface_dir = append(param_dir, "/interfaces");

        % Find the turnover points
        [ts, ds, Hs] = find_turnover_points(interface_dir, ...
            output_range, Delta, dt, transpose_coordinates, shift, ...
            curvature_time);

        % Save the turnover points as a text file
        combined = [ts; ds; Hs];
        writematrix(combined', sprintf("%s/turnover_points.csv", param_dir));
    end


%     cd(master_dir);
%     subDirs = dir;
% 
%     for k = 3 : length(subDirs)
%         currDir = subDirs(k).name
%     end
end