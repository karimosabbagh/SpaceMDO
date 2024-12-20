function [x, initial_guess] = generate_initial_guess()
    % This function generates a hardcoded initial guess for the problem definition (PB)
    % and saves it to a specified file. It also returns a 1x29 double array (x) 
    % containing the initial guesses.
    %
    % Inputs:
    % - file_name: The name of the file to save the initial guess (e.g., 'initial_guess.mat')
    %
    % Outputs:
    % - x: A 1x29 double array containing the initial guesses for each variable
    % - initial_guess: A structure containing the named initial guesses (saved to file)

    global R_earth R_mars start_date

    % Initialize the structure for storing the initial guess
    initial_guess = struct();

    % Hardcoded initial guesses based on the problem variables in PB.var
    initial_guess.delta_v_escape_e = 5.3;  % Example: mid-value within the bounds of [0, 80e3]
    initial_guess.r_p_e = 6.631;  %(R_earth + 400)/1e3;   % Example: some value between R_earth + 160e3 and R_earth + 2000e3
    %initial_guess.departure_date_e = start_date + 50; % Example: a date value within start_date to end_date
    %initial_guess.arrival_date_e = start_date + 200;  % Example: another date within the valid range
    initial_guess.V_SC_departure_e = 33.4;  % Example: value between 0 and 80e3
    initial_guess.V_SC_arrival_e = 19.7;    % Example: value between 0 and 80e3
    
    initial_guess.delta_v_capture_c = 9.8;  % Example: mid-value for delta_v_capture_c
    %initial_guess.m_prop_c = 1000;            % Example: some value for propellant mass
    initial_guess.r_p_c = 3.689; %(R_mars + 300)/1e3;    % Example: value between R_mars + 100e3 and 170 * R_mars
    initial_guess.e_c = 0.01;                 % Example: eccentricity value
    %initial_guess.departure_date_c = start_date + 60; % Example: departure date for capture
    %initial_guess.arrival_date_c = start_date + 180;  % Example: arrival date for capture
    initial_guess.V_SC_departure_c = 33.4;  % Example: speed value
    initial_guess.V_SC_arrival_c = 19.7;    % Example: speed value for arrival

    initial_guess.m_prop_s = 20000;            % Example propellant mass for spacecraft
    initial_guess.m_SC = 30000;      % Example structure mass
    initial_guess.Isp = 300;                 % Example ISP value
    initial_guess.delta_v_escape_s = 5.3;   % Example delta_v for spacecraft escape
    initial_guess.delta_v_capture_s = 9.8;  % Example delta_v for spacecraft capture
    %initial_guess.departure_date_s = start_date + 40; % Example: spacecraft departure date
    %initial_guess.arrival_date_s = start_date + 100;  % Example: spacecraft arrival date
    initial_guess.eta_FOV_s = deg2rad(15);           % Example: field-of-view angle for coverage
    initial_guess.IFOV_s = 10;                % Example: Instantaneous Field of View
    
    initial_guess.r_p_p = initial_guess.r_p_c;     % Example: planet coverage variable
    initial_guess.e_p = initial_guess.e_c;                % Example: eccentricity for planet coverage
    %initial_guess.T_orbit = 7000;            % Example: orbital period in seconds
    initial_guess.eta_center = deg2rad(0);  % Example: center angle for planet coverage
    initial_guess.eta_FOV = deg2rad(15);     % Example: field-of-view angle for coverage
    initial_guess.IFOV = initial_guess.IFOV_s;                  % Example: Instantaneous Field of View

    % Generate the 1x22 vector (x) based on the hardcoded initial guesses
    x = [
        initial_guess.delta_v_escape_e, initial_guess.r_p_e, initial_guess.V_SC_departure_e, initial_guess.V_SC_arrival_e, ...
        initial_guess.delta_v_capture_c, initial_guess.r_p_c, initial_guess.e_c, initial_guess.V_SC_departure_c, initial_guess.V_SC_arrival_c, ...
        initial_guess.m_prop_s, initial_guess.m_SC, initial_guess.Isp, initial_guess.delta_v_escape_s, initial_guess.delta_v_capture_s, initial_guess.eta_FOV_s, initial_guess.IFOV_s ...
        initial_guess.r_p_p, initial_guess.e_p, initial_guess.eta_center, initial_guess.eta_FOV, initial_guess.IFOV
    ];

    

    % Display a message to indicate the initial guess has been saved
    disp('Initial guess has been saved')
end
