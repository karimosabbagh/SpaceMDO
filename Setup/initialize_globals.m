function initialize_globals()
    % initialize_globals: Declare and initialize global variables

    global G M M_Sun R_earth R_mars earth_orbital_data mars_orbital_data date_range;

    % Define global constants
    G = 6.67430e-11;      % Gravitational constant (m^3/kg/s^2)
    M = 5.972e24;         % Mass of Earth (kg)
    M_Sun = 1.989e30;     % Mass of Sun (kg)
    R_earth = 6371e3;     % Radius of Earth (m)
    R_mars = 3389.5e3;    % Radius of Mars (m)

    start_date = datetime('2024-11-21', 'InputFormat', 'yyyy-MM-dd');
    end_date = datetime('2027-06-11', 'InputFormat', 'yyyy-MM-dd');
    date_step = days(1);
    date_range = start_date:date_step:end_date;
  
    % Load data from CSV files
    try
        earth_orbital_data = readtable('earth_result_table.csv');
        mars_orbital_data = readtable('mars_result_table.csv'); 

        % Display a confirmation message
        disp('Global variables and CSV data initialized successfully.');
    catch ME
        error('Error loading CSV files: %s', ME.message);
    end

end
