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

    % Calculate the objective function: total covered area over the elliptical orbit
    obj = -integrate_coverage(a, e, eta_center, eta_FOV_tilde);

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
    %A_total = integral(@(theta) instantaneous_coverage(a, e, theta, eta_center, eta_FOV_tilde), theta_min, theta_max);
    % Perform numerical integration
    try
        A_total = integral(@(theta) instantaneous_coverage(a, e, theta, eta_center, eta_FOV_tilde), theta_min, theta_max);
    catch ME
        fprintf("Error in integrate_coverage:\n  Message: %s\n", ME.message);
        A_total = NaN; % Return NaN to indicate failure
    end
    % Divide by 2*pi for average per anomaly (optional)
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

    % Instantaneous coverage area
    A = 2 * pi * R_mars^2 * (1 - cos(Lambda_FOV));
end
function Lambda = ground_range_angle(eta, R_planet, r_sat)
    % Solve for the ground range angle (Lambda) given eta
    % Inputs:
    %   eta        - Boresight angle (radians)
    %   R_planet   - Radius of the planet (m)
    %   r_sat      - Orbital radius (m)
    % Outputs:
    %   Lambda     - Ground range angle (radians)


    % Debugging step
    if ~isscalar(r_sat)
        disp('Error: r_sat is not scalar.');
        disp(size(r_sat)); % Show dimensions
        return; % Early exit for debugging
    end

    % Ensure r_sat is scalar
    assert(isscalar(r_sat), 'r_sat must be a scalar value');

    % Define the objective function for minimization
    ground_range_eq = @(Lambda) abs(tan(eta) - (R_planet * sin(Lambda)) / (r_sat - R_planet * cos(Lambda)));

    % Set bounds for Lambda [0, pi/2]
    lb = 0;           % Lower bound
    ub = pi / 2;      % Upper bound

    % Use eta as the initial guess
    Lambda0 = eta;

    % Solve using fmincon with bounds
    options = optimoptions('fmincon', 'Display', 'none'); % Suppress fmincon output
    Lambda = fmincon(ground_range_eq, Lambda0, [], [], [], [], lb, ub, [], options);
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
