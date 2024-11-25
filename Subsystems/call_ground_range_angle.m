% Wrapper function to call ground_range_angle from the command window
function Lambda = call_ground_range_angle(eta_deg, altitude)
    global R_mars;

    % Convert inputs
    eta = deg2rad(eta_deg);          % Convert degrees to radians
    R_planet = R_mars;               % Radius of Mars
    r_sat = R_planet + altitude;     % Orbital radius (m)

    % Call the ground_range_angle function
    Lambda = ground_range_angle(eta, R_planet, r_sat);

    % Display result (optional)
    fprintf('Ground range angle (Lambda): %.6f radians\n', Lambda);
end
