function initialize_globals()
    % initialize_globals: Declare and initialize global variables


    global G M_Earth M_Mars M_Sun R_earth R_mars Res_min earth_orbital_data mars_orbital_data start_date end_date departure_date arrival_date m_SC;
    
    % Define global constants
    G = 6.67430e-11 /1e9;      % Gravitational constant (m^3/kg/s^2)
    M_Earth = 5.972e24;   % Mass of Earth (kg)
    M_Mars = 0.6417e24;   % Mass of Mars (kg)
    M_Sun = 1.989e30;     % Mass of Sun (kg)
    R_earth = 6371;     % Radius of Earth (km)
    R_mars = 3389.5;    % Radius of Mars (km)

    % Mass of spacecraft

    m_SC = 3000; %kg

    % Planet Coverage Specific constant
    Res_min = 30; % m

    % Planetary Transfer date ranges
    start_date = datetime('2024-11-21', 'InputFormat', 'yyyy-MM-dd');
    start_date = datenum(start_date);
    end_date = datetime('2027-06-11', 'InputFormat', 'yyyy-MM-dd');
    end_date = datenum(end_date);

    % Fixed arrival and departure dates
    start_date = datetime('2024-11-21', 'InputFormat', 'yyyy-MM-dd');
    departure_date = datenum(start_date);
    departure_date = fix(departure_date);
    end_date = datetime('2025-05-20', 'InputFormat', 'yyyy-MM-dd');
    arrival_date = datenum(end_date);
    arrival_date = fix(arrival_date);
    
    % Load data from CSV files
    try
        earth_orbital_data = readtable('earth_result_table.csv');
        mars_orbital_data = readtable('mars_result_table.csv'); 


        % Convert the dates to datenums
        earth_orbital_data.DepartureDate = datetime(earth_orbital_data.DepartureDate, 'InputFormat', 'yyyy-MM-dd');
        earth_orbital_data.DateNum = datenum(earth_orbital_data.DepartureDate);
        mars_orbital_data.ArrivalDate = datetime(mars_orbital_data.ArrivalDate, 'InputFormat', 'yyyy-MM-dd');
        mars_orbital_data.DateNum = datenum(mars_orbital_data.ArrivalDate);



        % Display a confirmation message
        disp('Global variables and CSV data initialized successfully.');
    catch ME
        error('Error loading CSV files: %s', ME.message);
    end

    save globals M_Earth M_Mars M_Sun R_earth R_mars Res_min earth_orbital_data mars_orbital_data start_date end_date departure_date arrival_date m_SC;

end
