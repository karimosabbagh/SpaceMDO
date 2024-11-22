% define constants
M_Earth = 5.974e24;         % kg
M_Sun = 1.989e30;           % kg
r_Earth = 6378 * 1000;      % m
G = 6.6742e-11;             % m^3/kg*s^2

% constraint constants
r_p1_min = 80 * 1000;       % m
r_p1_max = 2000 * 1000;     % m
TOF_min = 128;              % days
TOF_max = 500;              % days

% inputs to vary in parametric study
% m_SC                       % kg (from Propellant and Structure Mass Subsystem)
% V_SC_arrival               % m/s (from Orbital Capture Subsystem)
% r_p1                       % m
% departure_date             % 'year-month-day'
% arrival_date               % 'year-month-day'


function [delta_v_escape, V_SC_departure] = orbital_escape_delta_v(m_SC, r_p1, V_SC_arrival, departure_date, arrival_date)

    % Extract Earth and Mars position and velocity for the exact departure and arrival dates
    departure_row = earth_result_table(earth_result_table.DepartureDate == departure_date, :);
    arrival_row = mars_result_table(mars_result_table.ArrivalDate == arrival_date, :);

    R_Earth_departure = [departure_row.Earth_Position_Magnitude] * 1000; % m
    V_Earth_departure = [departure_row.Earth_Velocity_Magnitude] * 1000; % m/s
    R_Mars_arrival = [arrival_row.Mars_Position_Magnitude] * 1000; % m



    % spacecraft departure velocity
    V_SC_departure = sqrt(VA_SC_arrival^2 + 2 * G * (M_Sun + m_SC) * ...
        ((1 / R_Earth_departure) - (1 / R_Mars_arrival)));
    
    % hyperbolic escape velocity
    V_infinity_D = V_SC_departure - V_Earth_departure;
    
    % calculate delta v
    delta_v_escape = sqrt(G * (M_Earth + m_SC) / (r_Earth + r_p1)) * ... 
        sqrt(2 + (V_infinity_D * sqrt(r_Earth + r_p1) / sqrt(G * (M_Earth + m_SC)))^2) - 1;

    
    % constraints

    % g1: Parking orbit constraint
    g1_lower = r_p1 - r_p1_min;                       % r_p1 >= r_p1_min
    g1_upper = r_p1_max - r_p1;                       % r_p1 <= r_p1_max

    % g2: Hyperbolic excess velocity constraint
    g2 = V_SC_departure - V_Earth;                    % V_D^(v) - V_Earth > 0

    % g3: Launch window constraint
    g3_lower = t_D - t_D_min;                         % t_D >= t_D_min
    g3_upper = t_D_max - t_D;                         % t_D <= t_D_max

    % g4: Time of flight constraint
    g4_tof_max = (t_A - t_D) - TOF_max;               % TOF_max >= t_A - t_D
    g4_tof_min = TOF_min - (t_A - t_D);               % t_A - t_D >= TOF_min

end