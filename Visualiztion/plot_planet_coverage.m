disp(['Integrated Coverage: ', num2str(A_total)]);
plot_coverage_2D(a, e, eta_center, eta_FOV_tilde)
Lambda = ground_range_angle(eta_center + eta_FOV_tilde / 2, R_mars, a);
plot_spherical_cap(Lambda);
plot_rotating_coverage3(a, e, eta_center, eta_FOV_tilde)

function plot_coverage_2D(a, e, eta_center, eta_FOV_tilde)
    global R_mars;

    % Define true anomaly range
    theta_range = linspace(0, 2*pi, 500); % Full orbit with 500 steps

    % Initialize arrays to store coverage limits
    lat_min = zeros(size(theta_range));
    lat_max = zeros(size(theta_range));

    % Loop over true anomaly
    for i = 1:length(theta_range)
        theta = theta_range(i);

        % Orbital radius at this true anomaly
        r_sat = a * (1 - e^2) / (1 + e * cos(theta));

        % Off-nadir angles
        eta_max = eta_center + eta_FOV_tilde / 2;
        eta_min = eta_center - eta_FOV_tilde / 2;

        % Ground range angles
        Lambda_max = ground_range_angle(eta_max, R_mars, r_sat);
        Lambda_min = ground_range_angle(eta_min, R_mars, r_sat);

        % Latitude limits of coverage
        lat_min(i) = -Lambda_min; % Southern edge of FoV
        lat_max(i) = Lambda_max;  % Northern edge of FoV
    end

    % Create a polar plot for coverage
    figure; % Open a new figure
    polaraxes; % Create polar axes
    hold on;

    % Plot northern and southern limits
    polarplot(theta_range, lat_max, 'b-', 'LineWidth', 1.5); % Northern edge
    polarplot(theta_range, lat_min, 'r-', 'LineWidth', 1.5); % Southern edge

    % Add labels and title
    title('Satellite Coverage in 2D Polar Plot');
    legend('Northern Coverage Limit', 'Southern Coverage Limit');

    hold off;
end


function plot_spherical_cap(Lambda)
    global R_mars;

    % Create a sphere
    [x, y, z] = sphere(100); % Smooth sphere
    x = x * R_mars;
    y = y * R_mars;
    z = z * R_mars;

    % Calculate the cap limit (z >= R*cos(Lambda))
    cap_limit = R_mars * cos(Lambda);
    mask = z >= cap_limit;

    % Create a mask for the cap
    cap_map = NaN(size(x));
    cap_map(mask) = 1; % Highlighted region

    % Plot the sphere with the cap
    figure;
    hold on;
    surface(x, y, z, 'FaceColor', 'texturemap', 'CData', cap_map, 'EdgeColor', 'none');
    colormap([0.8 0.8 0.8; 0 0.5 0]); % Gray for uncovered, green for covered
    colorbar;
    caxis([0 1]);
    title('Spherical Cap Coverage');
    axis equal;
    view(3);
    grid on;
    light('Position', [1 1 1], 'Style', 'infinite');
    lighting phong;
    hold off;
end

function plot_rotating_coverage(a, e, eta_center, eta_FOV_tilde)
    global R_mars;

    % Define true anomaly range
    theta_range = linspace(0, 2*pi, 100); % Full orbit with 100 steps
    % theta_range = [0, pi];

    % Angular radius of the field of view
    eta_max = eta_center + eta_FOV_tilde / 2;
    Lambda = ground_range_angle(eta_max, R_mars, a); % Fixed Lambda for circular orbit

    % Create a sphere
    [x, y, z] = sphere(100); % Smooth sphere
    x = x * R_mars;
    y = y * R_mars;
    z = z * R_mars;

    % Initialize coverage map
    coverage_map = NaN(size(x));

    % Loop through true anomaly to rotate the cap
    for theta = theta_range
        % Latitude and longitude of the sub-satellite point
        lat_center = eta_center; % Latitude offset by eta_center
        lon_center = theta;      % Longitude varies with orbital position

        % Rotate the cap
        for i = 1:size(x, 1)
            for j = 1:size(x, 2)
                % Calculate latitude and longitude for this grid point
                lat = asin(z(i, j) / R_mars); % Latitude
                lon = atan2(y(i, j), x(i, j)); % Longitude

                % Adjust longitude difference
                lon_diff = mod(lon - lon_center + pi, 2*pi) - pi; % Wrap to [-pi, pi]

                % Compute angular distance from the sub-satellite point
                angular_dist = acos(sin(lat_center) .* sin(lat) + ...
                    cos(lat_center) .* cos(lat) .* cos(lon_diff));

                % Check if this point is within the cap
                if angular_dist <= Lambda
                    coverage_map(i, j) = 1; % Mark as covered
                end
            end
        end
    end

    % Plot the sphere with the rotating coverage
    figure;
    hold on;
    surface(x, y, z, 'FaceColor', 'texturemap', 'CData', coverage_map, 'EdgeColor', 'none');
    colormap([0.8 0.8 0.8; 0 0.5 0]); % Gray for uncovered, green for covered
    colorbar;
    caxis([0 1]);
    title('Satellite Rotating Coverage in 3D');
    axis equal;
    view(3);
    grid on;
    light('Position', [1 1 1], 'Style', 'infinite');
    lighting phong;
    hold off;
end

function plot_rotating_coverage_2(a, e, eta_center, eta_FOV_tilde)
    global R_mars;

    % Define true anomaly range
    theta_range = linspace(0, 2*pi, 100); % Full orbit with 100 steps

    % Angular radius of the field of view
    eta_max = eta_center + eta_FOV_tilde / 2;
    eta_min = eta_center - eta_FOV_tilde / 2;

    % Ground range angles
    Lambda_max = ground_range_angle(eta_max, R_mars, a);
    Lambda_min = ground_range_angle(eta_min, R_mars, a);

    % Create a sphere
    [x, y, z] = sphere(100); % Smooth sphere
    x = x * R_mars;
    y = y * R_mars;
    z = z * R_mars;

    % Initialize coverage map
    coverage_map = NaN(size(x));

    % Loop through true anomaly to rotate the cap
    for theta = theta_range
        % Latitude and longitude of the sub-satellite point
        lat_center = eta_center; % Latitude offset by eta_center
        lon_center = theta;      % Longitude varies with orbital position

        % Rotate the cap
        for i = 1:size(x, 1)
            for j = 1:size(x, 2)
                % Calculate latitude and longitude for this grid point
                lat = asin(z(i, j) / R_mars); % Latitude
                lon = atan2(y(i, j), x(i, j)); % Longitude

                % Adjust longitude difference
                lon_diff = mod(lon - lon_center + pi, 2*pi) - pi; % Wrap to [-pi, pi]

                % Compute angular distance from the sub-satellite point
                angular_dist = acos(sin(lat_center) .* sin(lat) + ...
                    cos(lat_center) .* cos(lat) .* cos(lon_diff));

                % Check if this point is within the cap
                if angular_dist <= Lambda_max && angular_dist >= Lambda_min
                    coverage_map(i, j) = 1; % Mark as covered
                end
            end
        end
    end

    % Plot the sphere with the rotating coverage
    figure;
    hold on;
    surface(x, y, z, 'FaceColor', 'texturemap', 'CData', coverage_map, 'EdgeColor', 'none');
    colormap([0.8 0.8 0.8; 0 0.5 0]); % Gray for uncovered, green for covered
    colorbar;
    caxis([0 1]);
    title('Satellite Rotating Coverage with Lambda_min and Lambda_max');
    axis equal;
    view(3);
    grid on;
    light('Position', [1 1 1], 'Style', 'infinite');
    lighting phong;
    hold off;
end