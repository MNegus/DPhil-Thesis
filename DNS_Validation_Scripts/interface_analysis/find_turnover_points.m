function [ts, ds, Js] = find_turnover_points(interface_dir, output_range, ...
    Delta, dt, transpose_coordinates, shift, curvatureTime)

    %% Definitions
    tol = Delta; % Tolerance for extracting the interface
    ts = dt * (output_range - shift); 
    
    %% Magic numbers (later change to options)
    xMin = 2 * Delta;
    xMax = 0.5; % Maximum x value to search in
    
    %% Plot (optional)
%     close(figure(1));
%     figure(1);
    
    %% Boolean for switching to curvature
    curvature = false;

    %% Loop over timesteps
    for k = 1 : length(output_range)
        % Output timestep and time
        timestep = output_range(k)
        t = ts(k) 
        
        % Adjust minimum accordingly
        if timestep > 200
            xMin = 4 * Delta;
        end
        
        %% Read in interface file
        interface_filename = sprintf("%s/interface_%d.txt", ...
            interface_dir, timestep);
        [start_points, end_points] = ...
            read_interface_points(interface_filename, transpose_coordinates);

        %% Remove specific parts
        keepIdxs = (start_points(:, 1) < xMax) & (end_points(:, 1) < xMax) ...
            & (start_points(:, 1) > xMin) & (end_points(:, 1) > xMin);

        %% Extract interface
        [interface_start_points, interface_end_points] ...
            = extract_interface(start_points(keepIdxs, :), ...
            end_points(keepIdxs, :), tol);

        %% Find line segment centres
        seg_centres = 0.5 * (interface_start_points + interface_end_points);
        [xs, sort_idxs] = sort(seg_centres(:, 1));
        ys = seg_centres(sort_idxs, 2);
        
        %% Find the turnover point
        if (curvature == false)
            % If early, just take the minimum y value as the turover point
            [turnY, minIdx] = min(ys);
            turnX = xs(minIdx);
            % If turnY is greater than 1, then switch to curvature
            if (turnY >= 1)
                curvature = true;
            end
        else
            xMax = 0.25; % Maximum x value to search in
            % Else find where curvature is maximum
            d2y_dx2 = (ys(1 : end - 2) - 2 * ys(2 : end - 1) + ys(3 : end)) / (2 * Delta^2);
            dy_dx = (ys(3 : end) - ys(1 : end - 2)) / (2 * Delta^2);
            kappa = d2y_dx2 ./ (1 + dy_dx.^2).^(3/2);

            % Find maximum curvature
            [maxKappa, maxIdx] = max(kappa);
            if (isempty(xs(maxIdx)) || isempty(ys(maxIdx)))
                turnX = 0;
                turnY = 0;
            else
                turnX = xs(maxIdx);
                turnY = ys(maxIdx);
            end
        end
        
        % Save values
        turnX
        turnY
        ds(k) = turnY;
        Js(k) = turnX;
        
        % Plot
%         plot(xs, ys);
%         hold on;
%         scatter(turnX, turnY);
%         hold off;
%         drawnow;
        
    end
end