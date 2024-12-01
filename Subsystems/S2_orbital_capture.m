function [delta_v_capture, V_SC_arrival, S2_constraints] = S2_orbital_capture(e, delta_m_d, r_p2, V_SC_departure, departure_date, arrival_date)
% inputs to vary in parametric study
% m_SC                       % 3000 kg (from Propellant and Structure Mass Subsystem)
% delta_m_d                  % kg (from Propellant subsystem)
% V_SC_departure             % m/s (from Orbital Escape Subsystem)
% r_p2                       % m
% departure_date             % 'year-month-day'
% arrival_date               % 'year-month-day'

% Add subsystem paths
currentFilePath = fileparts(mfilename('fullpath'));
subsytems = fullfile(currentFilePath, '..', 'Setup');
addpath(subsytems);

 global G M_Mars M_Sun earth_orbital_data mars_orbital_data;

    % departure_date = datetime(departure_date, "ConvertFrom", "datenum", "Format", 'yyyy-MM-dd');
    % arrival_date = datetime(arrival_date, "ConvertFrom", "datenum", "Format", 'yyyy-MM-dd');
    
    departure_date = fix(departure_date);
    arrival_date = fix(arrival_date);

    % Extract Earth and Mars position and velocity for the exact departure and arrival dates
    departure_row = earth_orbital_data(earth_orbital_data.DateNum == departure_date, :);
    arrival_row = mars_orbital_data(mars_orbital_data.DateNum == arrival_date, :);

    R_Earth_departure = [departure_row.Earth_Position_Magnitude] * 1000; % m
    R_Mars_arrival = [arrival_row.Mars_Position_Magnitude] * 1000; % m
    V_Mars_arrival = [arrival_row.Mars_Velocity_Magnitude] * 1000; % m/s 
    m_SC = 3000; %kg

    % spacecraft departure velocity
    V_SC_arrival = sqrt(V_SC_departure^2 + 2 * G * (M_Sun + (m_SC - delta_m_d)) * ...
        ((1 / R_Mars_arrival) - (1 / R_Earth_departure)));
    
    % hyperbolic escape velocity
    V_infinity_A = V_SC_arrival - V_Mars_arrival;

    delta_v_capture = sqrt(V_infinity_A ^ 2 + 2 * G * (M_Mars + (m_SC - delta_m_d)) / (r_p2)) + ... 
        sqrt(G * (M_Mars + (m_SC - delta_m_d)) * (1 + e) / (r_p2));

    % Evaluate constraints
    S2_constraints = S2_evaluate_constraints(V_SC_arrival, V_Mars_arrival);
end 


function S2_constraints = S2_evaluate_constraints(V_SC_arrival, V_Mars_arrival)
                                                 
    % g1: Hyperbolic excess velocity constraint (V_SC_capture - V_Mars_arrival < 0)
    c1 =  V_Mars_arrival -  V_SC_arrival ;               

    % Combine constraints
    S2_constraints = c1;

end
