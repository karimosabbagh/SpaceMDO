% Subsystem function
function [obj, T_orbit, constraints] = S4_planet_coverage(r_p, e, eta_center, eta_FOV_tilde)
    % Planet Coverage Subsystem
    % Inputs:
    %   r_p          - Periapsis radius (m)
    %   e            - Eccentricity (dimensionless)
    %   eta_center   - Off-nadir angle (radians)
    %   eta_FOV_tilde - Field of View angle (radians)
    %
    % Outputs:
    %   obj          - Objective value (negative total coverage area for maximization)
    %   T_orbit      - Orbital period (s)
    %   constraints  - Array of constraint values [c1, c2]
    
    % Declare global variables
    global R_mars G M ;

    % Subsystem Specific Parameters
    global Res_min k;
    Res_min = 5;
    k = 1;

    % Calculate semi-major axis (a)
    a = r_p / (1 - e);

    % Calculate orbital period
    T_orbit = 2 * pi * sqrt(a^3 / (G * M));

    global i;
    i = 0;

    % Calculate total surface area of Mars  (assuming perfect sphere)
    A_mars = 4 * pi * R_mars ^2;
    % Calculate the objective function: total covered area over the elliptical orbit
    obj = -integrate_coverage(a, e, eta_center, eta_FOV_tilde) / A_mars;

    % Evaluate non-linear constraints
    constraints = evaluate_constraints(a, e, eta_center, eta_FOV_tilde);
end

function A_total = integrate_coverage(a, e, eta_center, eta_FOV_tilde)
    % Integrate coverage over the elliptical orbit
    global R_mars;

    % True anomaly range (0 to 2*pi)
    theta_min = 0;
    theta_max = 2 * pi;

    %% Perform numerical integration
    A_total = integral(@(theta) instantaneous_coverage(a, e, theta, eta_center, eta_FOV_tilde), theta_min, theta_max);
    % Perform numerical integration
    %try
    %    A_total = integral(@(theta) instantaneous_coverage(a, e, theta, eta_center, eta_FOV_tilde), theta_min, theta_max);
    
    %catch ME
    %    fprintf("Error in integrate_coverage:\n  Message: %s\n", ME.message);
    %    A_total = NaN; % Return NaN to indicate failure
    %end
    % Divide by 2*pi for average per anomaly (optional)
    disp(['Integrated Coverage: ', num2str(A_total)]);
    plot_coverage_2D(a, e, eta_center, eta_FOV_tilde)
    Lambda = ground_range_angle(eta_center + eta_FOV_tilde / 2, R_mars, a);
    plot_spherical_cap(Lambda);
    plot_rotating_coverage2(a, e, eta_center, eta_FOV_tilde)

end

function A = instantaneous_coverage(a, e, theta, eta_center, eta_FOV_tilde)
    % Calculate the instantaneous coverage area at a given true anomaly
    global R_mars;

    % Orbital radius at true anomaly
    r_sat = a * (1 - e^2) ./ (1 + e * cos(theta));

    disp('theta');
    disp(length(theta));
    % Calculate eta_max and eta_min
    eta_max = eta_center + eta_FOV_tilde / 2;
    eta_min = eta_center - eta_FOV_tilde / 2;

    Lambda_max = arrayfun(@(r) ground_range_angle(eta_max, R_mars, r), r_sat);
    Lambda_min = arrayfun(@(r) ground_range_angle(eta_min, R_mars, r), r_sat);

    % Effective field of view
    Lambda_FOV = Lambda_max - Lambda_min;


    % Save variables to a .mat file
    save('debug_lambda_values.mat', 'r_sat', 'Lambda_max', 'Lambda_min', 'Lambda_FOV');

    % Instantaneous coverage area
    % A = 2 * pi * R_mars^2 * (1 - cos(Lambda_FOV));
    S = R_mars * Lambda_FOV;
    A = S * R_mars; %   % Rate of change of area with respect to theta
    % dA_dtheta = S * R_mars; % Arc length times the band thickness factor
end

function Lambda = ground_range_angle(eta, R_planet, r_sat)
    % Calculate the ground range angle (Lambda) given eta
    % Inputs:
    %   eta        - Off-nadir angle (radians, can be negative)
    %   R_planet   - Radius of the planet (m)
    %   r_sat      - Orbital radius (m)
    % Outputs:
    %   Lambda     - Ground range angle (radians)

    % Ensure inputs are scalar
    assert(isscalar(r_sat), 'r_sat must be a scalar value');
    assert(isscalar(R_planet), 'R_planet must be a scalar value');
    assert(isscalar(eta), 'eta must be a scalar value');

    % Calculate gamma
    sin_gamma = (r_sat / R_planet) * sin(eta);
    gamma = asin(max(min(sin_gamma, 1), -1)); % Clamp value to [-1, 1]

    % Calculate slant range rho
    rho = R_planet * cos(gamma) + r_sat * cos(eta);

    % Calculate the ground range angle Lambda
    sin_Lambda = (rho / R_planet) * sin(eta);
    Lambda = asin(max(min(sin_Lambda, 1), -1)); % Clamp value to [-1, 1]

   
end


function constraints = evaluate_constraints(a, e, eta_center, eta_FOV_tilde)
    % Evaluate the non-linear constraints
    global R_mars Res_min k;

    % Orbital radius at periapsis
    r_sat = a * (1 - e^2) / (1 + e * cos(0)); % r_sat for theta = 0

    % Constraint 1: Field of view does not extend beyond the planet's surface
    eta_max = eta_center + eta_FOV_tilde / 2;
    c1 = sin(eta_max) - R_mars / r_sat;

    % Constraint 2: Minimum resolution constraint
    Lambda_max = ground_range_angle(eta_max, R_mars, r_sat);

    rho = R_mars / sin(Lambda_max); % Slant range
    resolution = rho * k * cos(eta_center);
    c2 = Res_min - resolution;

    % Combine constraints
    constraints = [c1, c2];
end
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
        % Longitude of the sub-satellite point
        lon_center = theta;

        % Rotate the cap
        for i = 1:size(x, 1)
            for j = 1:size(x, 2)
                % Calculate latitude and longitude for this grid point
                lat = asin(z(i, j) / R_mars); % Latitude
                lon = atan2(y(i, j), x(i, j)); % Longitude

                % Adjust longitude difference
                lon_diff = mod(lon - lon_center + pi, 2*pi) - pi; % Wrap to [-pi, pi]

                % Check if this point is within the cap
                if lat >= -Lambda && lat <= Lambda && abs(lon_diff) <= Lambda
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

function plot_rotating_coverage2(a, e, eta_center, eta_FOV_tilde)
    global R_mars;

    % Define true anomaly range
    theta_range = linspace(0, 2*pi, 100); % Full orbit with 100 steps

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
