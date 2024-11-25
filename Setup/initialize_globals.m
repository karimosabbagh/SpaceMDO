function initialize_globals()
    % initialize_globals: Declare and initialize global variables

    global G M_Earth M_Mars M_Sun R_earth R_mars earth_orbital_data mars_orbital_data;

    % Define global constants
    G = 6.67430e-11;      % Gravitational constant (m^3/kg/s^2)
    M_Earth = 5.972e24;         % Mass of Earth (kg)
    M_Mars = 0.6417e24;
    M_Sun = 1.989e30;     % Mass of Sun (kg)
    R_earth = 6371e3;     % Radius of Earth (m)
    R_mars = 3389.5e3;    % Radius of Mars (m)
  
    % Load data from CSV files
    try
        earth_orbital_data = readmatrix('earth_result_table.csv');
        mars_orbital_data = readmatrix('mars_result_table.csv'); 

        % Display a confirmation message
        disp('Global variables and CSV data initialized successfully.');
    catch ME
        error('Error loading CSV files: %s', ME.message);
    end

end
