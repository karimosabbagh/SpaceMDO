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


function [delta_v_capture, V_SC_arrival] = orbital_capture_delta_v(e, m_SC,delta_m_d, r_p2, V_SC_departure, departure_date, arrival_date)

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
    e = ((2 * G * (M_mars + (m_SC - delta_m_d) / V_infinity_A^2) - r_p2) / ...
        (2 * G * (M_mars + (m_SC - delta_m_d) / V_infinity_A^2) + r_p2);
    
    % calculate delta v
    delta_v_capture = sqrt(V_infinity_A ^ 2 + 2 * G * (M_mars + (m_SC - delta_m_d)) / (r_p2)) + ... 
        sqrt(G * (M_mars + (m_SC - delta_m_d)) * (1 + e) / (r_p2));


    % constraints

    % g1: Parking orbit constraint
    g1_lower = r_p2 - r_p2_min;                       % r_p2 >= r_p2_min
    g1_upper = r_p2_max - r_p2;                       % r_p2 <= r_p2_max

    % g2: Hyperbolic excess velocity constraint
    g2 = V_SC_arrival - V_Mars_arrival;                    % V_A^(v) - V_Mars_arrival > 0

    % g3: Launch window constraint
    % implicit in earth_result_table as the possible date ranges are already provided

    % g4: Time of flight constraint
    g4_tof_max = (t_A - t_D) - TOF_max;               % TOF_max >= t_A - t_D
    g4_tof_min = TOF_min - (t_A - t_D);               % t_A - t_D >= TOF_min

    %g5: eccentricity constraint                      % 0 < e < 1
    g5 

end
