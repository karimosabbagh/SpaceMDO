% subsystem specific constraint parameters
r_p1_min = 160 * 1000;      % (m)
r_p1_max = 2000 * 1000;     % (m)
TOF_min = 128;              % (days)
TOF_max = 500;              % (days)


function [obj, V_SC_departure, constraints] = orbital_escape_delta_v(m_SC, r_p1, V_SC_arrival, departure_date, arrival_date)
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
    %   constraints       - Array of constraint values [c1, c2]

    global R_earth G M M_Sun earth_orbital_data mars_orbital_data;

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
    delta_v_escape = sqrt(G * (M + m_SC) / (R_earth + r_p1)) * ... 
        sqrt(2 + (V_infinity_D * sqrt(R_earth + r_p1) / sqrt(G * (M + m_SC)))^2) - 1;

    obj = delta_v_escape;
    
    % Evaluate constraints
    constraints = evaluate_constraints(r_p1, r_p1_min, r_p1_max, V_SC_departure, V_Earth_departure, arrival_date, departure_date);

end


function constraints = evaluate_constraints(r_p1, r_p1_min, r_p1_max, V_SC_departure, V_Earth_departure, arrival_date, departure_date)

    % g1: Parking orbit constraint
    c1 = r_p1 - r_p1_min;                                         % r_p1 >= r_p1_min
    c2 = r_p1_max - r_p1;                                         % r_p1 <= r_p1_max

    % g2: Hyperbolic excess velocity constraint
    c3 = V_SC_departure - V_Earth_departure;                      % V_D^(v) - V_Earth > 0

    % g3: Launch window constraint
    % implicit in earth_orbital_data as the possible date ranges are already provided

    % g4: Time of flight constraint
    arrival_date_datetime = datetime(arrival_date, 'InputFormat', 'yyyy-MM-dd');
    departure_date_datetime = datetime(departure_date, 'InputFormat', 'yyyy-MM-dd');
    date_diff = arrival_date_datetime - departure_date_datetime;

    c4 = date_diff - TOF_max;               % TOF_max >= t_A - t_D
    c5 = TOF_min - date_diff;               % t_A - t_D >= TOF_min

    % Combine constraints
    constraints = [c1, c2, c3, c4, c5];

end