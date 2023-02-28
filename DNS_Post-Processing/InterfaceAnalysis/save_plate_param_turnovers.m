function save_plate_param_turnovers(master_dir)
    close all;

    %% Options
    max_level = 15;
    output_range = 0 : 700;
    box_width = 24;
    Delta = box_width ./ 2^max_level;
    dt = 1e-3;
    transpose_coordinates = false;
    shift = 0;
    curvature_time = 0.45;

    %% Loop over varying cases
%     params = ["ALPHA", "BETA", "GAMMA"];
%     params = "GAMMA";
    params = "BETA_ALTERNATIVE";

    for param = params
%         parent_dir = sprintf("%s/%s_varying/", master_dir, param);
%         parent_dir = "/scratch/negus/MembraneRubberRuns/ALPHA_1.1_GAMMA_668";

        parent_dir = "/scratch/negus/MembraneStationaryCase"

        %% Loop over subdirectories
%         dirStruct = dir(parent_dir);

%         for k = 3 : length(dirStruct)
%             param_dir = append(parent_dir, dirStruct(k).name)
            param_dir = parent_dir;

            % Interface directory
            interface_dir = append(param_dir, "/interfaces");

            % Find the turnover points
            [ts, ds, Hs] = find_turnover_points(interface_dir, ...
                output_range, Delta, dt, transpose_coordinates, shift, ...
                curvature_time);

            % Save the turnover points as a text file
            combined = [ts; ds; Hs];
            writematrix(combined', sprintf("%s/turnover_points.csv", param_dir));
%         end
    end


%     cd(master_dir);
%     subDirs = dir;
% 
%     for k = 3 : length(subDirs)
%         currDir = subDirs(k).name
%     end
end