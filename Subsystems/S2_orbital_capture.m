% define constants
M_Earth = 5.974e24;         % kg
M_Sun = 1.989e30;           % kg
M_mars = 6.417e23;          % kg
r_Earth = 6378 * 1000;      % m
r_mars = 3396 * 1000;       % m
G = 6.6742e-11;             % m^3/kg*s^2

% constraint constants
r_p2_min = 150 * 1000;       % m
r_p2_max = 1000 * 1000;     % m
TOF_min = 128;              % days
TOF_max = 500;              % days

% inputs to vary in parametric study
% m_SC                       % kg (from Propellant and Structure Mass Subsystem)
% delta_m_d                  % kg (from Propellant subsystem)
% V_SC_departure             % m/s (from Orbital Escape Subsystem)
% r_p2                       % m
% departure_date             % 'year-month-day'
% arrival_date               % 'year-month-day'
function [delta_v_capture, V_SC_arrival, S2_constraint] = orbital_capture_delta_v(e, m_SC, delta_m_d, r_p2, V_SC_departure, departure_date, arrival_date)

 global R_earth G M_Earth M_Sun earth_orbital_data mars_orbital_data;

    % Extract Earth and Mars position and velocity for the exact departure and arrival dates
    departure_row = earth_result_table(earth_result_table.DepartureDate == departure_date, :);
    arrival_row = mars_result_table(mars_result_table.ArrivalDate == arrival_date, :);

    R_Earth_departure = [departure_row.Earth_Position_Magnitude] * 1000; % m
    R_Mars_arrival = [arrival_row.Mars_Position_Magnitude] * 1000; % m
    V_Mars_arrival = [arrival_row.Mars_Velocity_Magnitude] * 1000; % m/s


    % spacecraft departure velocity
    V_SC_arrival = sqrt(V_SC_departure^2 + 2 * G * (M_Sun + (m_SC - delta_m_d)) * ...
        ((1 / R_Mars_arrival) - (1 / R_Earth_departure)));
    
    % hyperbolic escape velocity
    V_infinity_A = V_SC_arrival - V_Mars_arrival;

        % calculate eccentricity
    e = ((2 * G * (M_mars + (m_SC - delta_m_d)) / V_infinity_A^2) - r_p2) / ...
        ((2 * G * (M_mars + (m_SC - delta_m_d)) / V_infinity_A^2) + r_p2);
    
    % calculate delta v
    delta_v_capture = sqrt(V_infinity_A ^ 2 + 2 * G * (M_mars + (m_SC - delta_m_d)) / (r_p2)) + ... 
        sqrt(G * (M_mars + (m_SC - delta_m_d)) * (1 + e) / (r_p2));

    % Evaluate constraints
    S2_constraints = S2_evaluate_constraints(r_p2, e, V_SC_arrival, V_Mars_arrival, ...
        departure_date, arrival_date);
end 

function tof = determine_tof(departure_date, arrival_date)
% Compute the Julian day number at 0 UT for departure and arrival dates, and the time of flight
    
    % Extract year, month, and day
    departure_year = year(departure_date);
    departure_month = month(departure_date);
    departure_day = day(departure_dt);

    arrival_year = year(arrival_date);
    arrival_month = month(arrival_date);
    arrival_day = day(arrival_date);
    
    % Compute J0 values
    jd1 = 367*departure_year - fix(7*(departure_year + fix((departure_month + 9)/12))/4) + ...
        fix(275*departure_month/9) + departure_day + 1721013.5;
    jd2 = 367*arrival_year - fix(7*(arrival_year + fix((arrival_month + 9)/12))/4) + ...
        fix(275*arrival_month/9) + arrival_day + 1721013.5;

    % Determine TOF
    tof = jd2 - jd1;
end

function S2_constraints = S2_evaluate_constraints(r_p2, V_SC_arrival, V_Mars_arrival, ...
    departure_date, arrival_date)

  % g1: Parking orbit constraint (r_p2_min <= r_p2 <= r_p2_max)
    c1 = r_p2;                                                     

    % g2: Hyperbolic excess velocity constraint (V_D^(v) - V_Earth > 0)
    c2 = V_SC_capture - V_Mars_arrival;               

    % g3: Launch window constraint (launch window start <= departure_date <= launch window end)
    c3 = departure_date;
    
    % g4: Time of flight constraint; (tof_min <= t_A - t_D <= tof_max)
    tof = determine_tof(departure_date, arrival_date);
    c4 = tof;              

    % g5: Eccentricity
    c5 = e;

    % Combine constraints
    S1_constraints = [c1, c2, c3, c4, c5];

end