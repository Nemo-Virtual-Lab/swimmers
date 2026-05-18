function [pointsInsideIntersection, pointsOutsideIntersection] = pointsInsideIntersection(convexHulls, numPoints, range)
    % Generate random points in the specified range
    pointsToCheck = bsxfun(@times, rand(numPoints, 3), range(:, 2)' - range(:, 1)') + range(:, 1)';
    

    % Initialize logical array to track points inside each convex hull
    inside = false(numPoints, numel(convexHulls));
    
    % Check points against each convex hull
    for i = 1:numel(convexHulls)
        % Determine which points are inside the current convex hull
        inside(:, i) = inpolyhedron(convexHulls{i}, pointsToCheck);
    end

    % Determine which points are inside at least 2 convex hulls
    insideAtLeastTwo = sum(inside, 2) >= 2;
    
    % Separate points into inside and outside groups of the intersection
    pointsInsideIntersection = pointsToCheck(insideAtLeastTwo, :);
    pointsOutsideIntersection = pointsToCheck(~insideAtLeastTwo, :);
end
