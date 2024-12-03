function [delta_v_capture, V_SC_arrival, V_infinity_A] = S2_orbital_capture(e, m_prop, r_p2, V_SC_departure)
% inputs to vary in parametric study
% m_SC                       % 3000 kg (from Propellant and Structure Mass Subsystem)
% m_prop                  % kg (from Propellant subsystem)
% V_SC_departure             % m/s (from Orbital Escape Subsystem)
% r_p2                       % m
% departure_date             % 'year-month-day'
% arrival_date               % 'year-month-day'

% Add subsystem paths
currentFilePath = fileparts(mfilename('fullpath'));
subsytems = fullfile(currentFilePath, '..', 'Setup');
addpath(subsytems);

 global G M_Mars M_Sun earth_orbital_data mars_orbital_data m_SC departure_date arrival_date;

 % Extract Earth and Mars position and velocity for the exact departure and arrival dates
 departure_row = earth_orbital_data(earth_orbital_data.DateNum == departure_date, :);
 arrival_row = mars_orbital_data(mars_orbital_data.DateNum == arrival_date, :);

 R_Earth_departure = [departure_row.Earth_Position_Magnitude] * 1000; % m
 R_Mars_arrival = [arrival_row.Mars_Position_Magnitude] * 1000; % m
 V_Mars_arrival = [arrival_row.Mars_Velocity_Magnitude] * 1000; % m/s

 % spacecraft departure velocity
 V_SC_arrival = sqrt(V_SC_departure^2 + 2 * G * (M_Sun + (m_SC - m_prop)) * ...
     ((1 / R_Mars_arrival) - (1 / R_Earth_departure)));
    
 % hyperbolic escape velocity
 V_infinity_A = V_Mars_arrival -  V_SC_arrival;

 delta_v_capture = sqrt(V_infinity_A ^ 2 + 2 * G * (M_Mars + (m_SC - m_prop)) / (r_p2)) + ... 
     sqrt(G * (M_Mars + (m_SC - m_prop)) * (1 + e) / (r_p2));

end 
