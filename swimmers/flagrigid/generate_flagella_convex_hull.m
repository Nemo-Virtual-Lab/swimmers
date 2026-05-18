function [vtx, elt] = generate_flagella_convex_hull(N, L, rad, step_length, rad_section, ke_type, nb_rot)


    if strcmp(ke_type, 'schum')
        k_E=0.333*2*pi/step_length;% Schum
    elseif strcmp(ke_type, 'pt')
        k_E=2*pi/step_length; % Phan-Thien
    end

    % Generate points along the flagella
    S = linspace(0, L, N);
    x = S;
    y = rad*(1 - exp(-k_E^2 * S.^2));
    z = rad*(1 - exp(-k_E^2 * S.^2)) + rad_section;
    points = [x', y', z'];

    % Rotate the points around the axis
    theta = 0;
    L_points = {points};
    for i = 1:nb_rot-1
        % Define rotation parameters
        theta = theta + 2*pi / nb_rot; % Angle of rotation
        axis = [1, 0, 0];  % Rotation axis, you can change it as needed

        % Rotate the points
        rotation_matrix = [cos(theta) + axis(1)^2 * (1 - cos(theta)),...
                           axis(1) * axis(2) * (1 - cos(theta)) - axis(3) * sin(theta),...
                           axis(1) * axis(3) * (1 - cos(theta)) + axis(2) * sin(theta);...
                           axis(2) * axis(1) * (1 - cos(theta)) + axis(3) * sin(theta),...
                           cos(theta) + axis(2)^2 * (1 - cos(theta)),...
                           axis(2) * axis(3) * (1 - cos(theta)) - axis(1) * sin(theta);...
                           axis(3) * axis(1) * (1 - cos(theta)) - axis(2) * sin(theta),...
                           axis(3) * axis(2) * (1 - cos(theta)) + axis(1) * sin(theta),...
                           cos(theta) + axis(3)^2 * (1 - cos(theta))];

        rotated_points = points * rotation_matrix;
        L_points{end+1} = rotated_points;
    end

    % Compute convex hull
    N_points = vertcat(L_points{:});
    k = convhull(N_points(:,1), N_points(:,2), N_points(:,3));

    % Extract vertices and faces of convex hull
    vtx = N_points;
    elt = k;
end
