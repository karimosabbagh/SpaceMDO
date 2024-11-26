function [delta_v_escape, y , S1_constraints] = ...
    orbital_escape_delta_v(m_SC, r_p1, V_SC_arrival, departure_date, arrival_date)
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
    %   constraints       - Array of constraint values [c1, c2, c3, c4]

    global G M_Earth M_Sun earth_orbital_data mars_orbital_data;

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
        sqrt(2 + (V_infinity_D * sqrt( r_p1) / sqrt(G * (M_Earth + m_SC)))^2) - 1;

    tof = determine_tof(departure_date, arrival_date);

    y = [V_SC_departure tof];

    % Evaluate constraints
    S1_constraints = S1_evaluate_constraints(r_p1, V_SC_departure, V_Earth_departure, ...
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


function S1_constraints = S1_evaluate_constraints(V_SC_departure, V_Earth_departure)

   % g1: Hyperbolic excess velocity constraint (V_D^(v) - V_Earth > 0)
    c1 = V_SC_departure - V_Earth_departure;               
              
    % Combine constraints
    S1_constraints = c1;

end