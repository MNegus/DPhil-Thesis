function [interface_start_points, interface_end_points] ...
    = extract_interface(start_points, end_points, tol)
%extract_interface.m Extracts the bulk droplet interface from the
%collection of all of the VOF points.
%   We assume the coordinates are setup such that the vertical coordinate
%   is along x and the radial/horizontal coordinate along y.

    %% Parameters
%     tol = 0.75 * 6 / 2^11; % Minimum distance for a cell to be seen as adjacent

    %% Collected arrays
%     % Collection of all of the segments, with unique numbering
    
%     segs = zeros(noSegs, 5);
%     segs(:, 1) = 1 : noSegs;
%     segs(:, 2 : 3) = start_points;
%     segs(:, 4 : 5) = end_points;

    % Collection of all of the points, with unique line segment numbering, and
    % numbering indicators whether the point is start or end, with 1 for start
    % and -1 for end
    noSegs = length(start_points);
    allPoints = zeros(2 * noSegs, 4);

    allPoints(1 : noSegs, 1) = 1 : noSegs;
    allPoints(1 : noSegs, 2 : 3) = start_points;
    allPoints(1 : noSegs, 4) = ones(noSegs, 1);

    allPoints(noSegs + 1 : end, 1) = 1 : noSegs;
    allPoints(noSegs + 1 : end, 2 : 3) = end_points;
    allPoints(noSegs + 1 : end, 4) = -ones(noSegs, 1);


    %% Function to find the opposite point on a line segment
    function  point = oppositePoint(allPoints, segNo, m)
       % Find the point that has the same segNo but opposite m
       point ...
           = allPoints((allPoints(:, 1) == segNo) ...
            & (allPoints(:, 4) == -m), 2 : 3);
    end

    %% Find the first segment, which is the one at the top of the droplet
    % This segment will either start or end with y = 0.

    % Find maximum x value
    [xMax, maxIdx] = max(allPoints(:, 2));

    % Save the first start point
    maxSegNo = allPoints(maxIdx, 1);
    initStart = allPoints(maxIdx, 2 : 3);
    interface_start_points = [initStart];
    startm = allPoints(maxIdx, 4);

    % Save the first end point
    initEnd = oppositePoint(allPoints, maxSegNo, startm);
    interface_end_points = initEnd;
    
    % Remove points from the array
    allPoints(allPoints(:, 1) == maxSegNo, :) = [];
    
    %% Loop over the points until we've found all the connected points
    currPoint = initEnd;
    currm = -startm;
    closestDist = 0;
    loopCounter = 0;
    
    while (closestDist < tol) && (loopCounter < noSegs - 1)
        
        %% Distances from any point to the current point
        dists = sqrt((allPoints(:, 2) - currPoint(1)).^2 ...
            + (allPoints(:, 3) - currPoint(2)).^2);
        
        %% Find minimum distance
        [closestDist, minIdx] = min(dists);
        
        if (closestDist < tol) 
           % Find the segment from allPoints
           segNo = allPoints(minIdx, 1);
           startPoint = allPoints(minIdx, 2 : 3);
           m = allPoints(minIdx, 4);
           endPoint = oppositePoint(allPoints, segNo, m);
           
           % Save the segment
           interface_start_points(end + 1, :) = startPoint;
           interface_end_points(end + 1, :) = endPoint;
           
           % Remove segment from allPoints
           allPoints(allPoints(:, 1) == segNo, :) = [];
           
           % Update currPoint
           currPoint = endPoint;
        end
        
        % Increment loop counter
        loopCounter = loopCounter + 1;
    end
    
    
    %% Plot

%     scatter(allPoints(:, 2), allPoints(:, 3));
%     scatter(interface_start_points(:, 1), interface_start_points(:, 2));

end