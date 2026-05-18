function [is_intersection, intersection_point] = triangleHalfLineIntersection(...
    half_line_origin, half_line_direction, vertex1, vertex2, vertex3)

    triangle_vertices = [vertex1; vertex2; vertex3];
    
    % Compute normal vector of the triangle
    v1 = vertex2 - vertex1;
    v2 = vertex3 - vertex1;
    normal_vector = cross(v1,v2);
    
    % Check if the half-line intersects the plane of the triangle
    if dot(normal_vector, half_line_direction) == 0
        is_intersection = false;
        intersection_point = [];
        return;
    end
    
    % Compute the intersection point in the plane of the triangle
    t = dot(normal_vector, (vertex1 - half_line_origin)) / dot(normal_vector, half_line_direction);
    intersection_point_plane = half_line_origin + t * half_line_direction;
    
    % Check if the intersection point is inside the triangle
    if isInsideTriangle(intersection_point_plane, triangle_vertices)
        is_intersection = true;
        intersection_point = [intersection_point_plane]; % Add a placeholder for the z-coordinate
    else
        is_intersection = false;
        intersection_point = [];
    end
end