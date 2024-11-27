function [delta_v_escape, V_SC_departure, S1_constraints] = ...
    S1_orbital_escape(m_SC, r_p1, V_SC_arrival, departure_date, arrival_date)
    %
    % Determine the delta_v to escape Earth's gravitational sphere of influence and 
    % velocity of spacecraft at departure
    %
    % Inputs:
    %   m_SC              - spacecraft mass (kg)
    %   r_p1              - radius of parking orbit (m)
    %   V_SC_arrival      - velocity of spacecraft at Mars arrival SOI (m/s)
    %   departure_date    - date of departure (year-month-day)
    %   arrival_date      - date of arrival (year-month-day)
    %
    % Outputs:
    %   obj               - objective value (delta_v to escape Earth SOI (m/s))
    %   V_SC_departure    - velocity of spacecraft at Earth SOI (m/s)
    %   constraints       - Array of constraint values [c1, c2, c3]

    % Add subsystem paths
    currentFilePath = fileparts(mfilename('fullpath'));
    setup = fullfile(currentFilePath, '..', 'Setup');
    addpath(setup);
    

    global G M_Earth M_Sun earth_orbital_data mars_orbital_data;

    departure_date = datetime(departure_date, "ConvertFrom", "posixtime", "Format", 'yyyy-MM-dd');
    arrival_date = datetime(arrival_date, "ConvertFrom", "posixtime", "Format", 'yyyy-MM-dd');

    % Extract Earth and Mars position and velocity for the exact departure and arrival dates
    departure_row = earth_orbital_data(earth_orbital_data.DepartureDate == departure_date, :);
    arrival_row = mars_orbital_data(mars_orbital_data.ArrivalDate == arrival_date, :);
    R_Earth_departure = [departure_row.Earth_Position_Magnitude] * 1000; % m
    V_Earth_departure = [departure_row.Earth_Velocity_Magnitude] * 1000; % m/s
    R_Mars_arrival = [arrival_row.Mars_Position_Magnitude] * 1000; % m

    % spacecraft departure velocity
    V_SC_departure = sqrt(V_SC_arrival^2 + 2 * G * (M_Sun + m_SC) * ...
        ((1 / R_Earth_departure) - (1 / R_Mars_arrival)));
    
    % hyperbolic escape velocity
    V_infinity_D = V_SC_departure - V_Earth_departure;
    
    % calculate delta v
    delta_v_escape = sqrt(G * (M_Earth + m_SC) / (r_p1)) * ... 
        sqrt(2 + (V_infinity_D * sqrt(r_p1) / sqrt(G * (M_Earth + m_SC)))^2) - 1;

    % calculate time of flight
    tof = determine_tof(departure_date, arrival_date);

    % Evaluate constraints
    S1_constraints = S1_evaluate_constraints(V_SC_departure, V_Earth_departure, tof);

end


function S1_constraints = S1_evaluate_constraints(V_SC_departure, V_Earth_departure, tof)

   % g1: Hyperbolic excess velocity constraint (V_D^(v) - V_Earth > 0)
    c1 = V_Earth_departure - V_SC_departure;   

   % g2: Time of flight
    c2 = 128 - tof;
    c3 = tof - 500;
              
    % Combine constraints
    S1_constraints = [c1, c2, c3];

end