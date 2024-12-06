% Subsystem function
function [obj, S4_constraints] = S4_planet_coverage(r_p, e, eta_center, eta_FOV_tilde, IFOV)
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
    global R_mars G M_Mars;

    % Subsystem Specific Parameters
    global Res_min;

    IFOV = IFOV/1000;

    % Calculate semi-major axis (a)
    a = r_p * 1e3 / (1 - e); % Convert to km

    % Calculate orbital period
    T_orbit = 2 * pi * sqrt(a^3 / (G * M_Mars));

    % Calculate total surface area of Mars  (assuming perfect sphere)
    A_mars = 4 * pi * R_mars ^2;
    % Calculate the objective function: total covered area over the elliptical orbit
    obj = -integrate_coverage(a, e, eta_center, eta_FOV_tilde) / A_mars;

    % Evaluate non-linear constraints
    S4_constraints = evaluate_constraints(a, e, eta_center, eta_FOV_tilde, IFOV);
end

function A_total = integrate_coverage(a, e, eta_center, eta_FOV_tilde)
    % Integrate coverage over the elliptical orbit
    global R_mars;

    % True anomaly range (0 to 2*pi)
    theta_min = 0;
    theta_max = 2 * pi;

    %% Perform numerical integration
    % A_total = integral(@(theta) instantaneous_coverage(a, e, theta, eta_center, eta_FOV_tilde), theta_min, theta_max);
    % Number of points for numerical integration (higher = more accurate)
    N_points = 500;
    
    % Discretize the range of theta
    theta_values = linspace(theta_min, theta_max, N_points);
    
    % Evaluate the integrand at each point
    integrand_values = arrayfun(@(theta) instantaneous_coverage(a, e, theta, eta_center, eta_FOV_tilde), theta_values);

    % Apply the trapezoidal rule for numerical integration
    A_total = trapz(theta_values, integrand_values);

end

function A = instantaneous_coverage(a, e, theta, eta_center, eta_FOV_tilde)
    % Calculate the instantaneous coverage area at a given true anomaly
    global R_mars;

    % Orbital radius at true anomaly
    r_sat = a * (1 - e^2) ./ (1 + e * cos(theta));

    
    % Calculate eta_max and eta_min
    eta_max = eta_center + eta_FOV_tilde / 2;
    eta_min = eta_center - eta_FOV_tilde / 2;

    Lambda_max = arrayfun(@(r) ground_range_angle(eta_max, R_mars, r), r_sat);
    Lambda_min = arrayfun(@(r) ground_range_angle(eta_min, R_mars, r), r_sat);

    % Effective field of view
    Lambda_FOV = Lambda_max - Lambda_min;


    % Save variables to a .mat file
    % save('debug_lambda_values.mat', 'r_sat', 'Lambda_max', 'Lambda_min', 'Lambda_FOV');

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


function constraints = evaluate_constraints(a, e, eta_center, eta_FOV_tilde, IFOV)
    % Evaluate the non-linear constraints
    global R_mars Res_min;

    % Orbital radius at periapsis
    r_sat = a * (1 - e^2) / (1 + e * cos(pi)); % r_sat for theta = 180

    % Constraint 1: Field of view does not extend beyond the planet's surface
    eta_max = eta_center + eta_FOV_tilde / 2;
    c1 = sin(eta_max) - R_mars / r_sat;

    % Constraint 2: Minimum resolution constraint
    % Lambda_max = ground_range_angle(eta_max, R_mars, r_sat);
    % rho = R_mars / sin(Lambda_max); % Slant range
    % resolution = rho * k * cos(eta_center);
    % c2 = resolution - Res_min ;
    

    % Calculate gamma
    sin_gamma = (r_sat / R_mars) * sin(eta_max);
    
    gamma = pi - asin(max(min(sin_gamma, 1), -1)); % Clamp value to [-1, 1] get oblique angle
    
    % Calculate slant range rho
    rho = R_mars * cos(gamma) + r_sat * cos(eta_max);
    
    c2 = rho * IFOV - Res_min; %rho is in km, IFOV in mrad

    if isnan(c1)
        disp('c1 is NaN')
    end

    if isnan(c2)
        disp('c2 is NaN')
    end

    % Combine constraints
    constraints = [c1, c2];
end

