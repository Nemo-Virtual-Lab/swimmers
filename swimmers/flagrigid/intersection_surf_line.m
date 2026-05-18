function [intersection_points] = intersection_surf_line(x_s, mesh, alpha, beta)
% Define the half-line in xOy plane passing through x_s and pointing along alpha
half_line_origin = x_s'; % Take only x and y coordinates
half_line_direction = [cos(alpha)*cos(beta), sin(beta)*cos(alpha), sin(alpha)];
% Check for intersection with each triangle in the mesh
intersection_points = [];
for i = 1:size(mesh.elt, 1)
    vertex1 = mesh.vtx(mesh.elt(i,1),:);
    vertex2 = mesh.vtx(mesh.elt(i,2),:);
    vertex3 = mesh.vtx(mesh.elt(i,3),:);
    % Check if the triangle intersects with the half-line
    [is_intersection, intersection_point] = triangleHalfLineIntersection(...
        half_line_origin, half_line_direction, vertex1, vertex2, vertex3);
    % If there is an intersection, store the point
    if is_intersection
        intersection_points = [intersection_points; intersection_point];
    end
end
