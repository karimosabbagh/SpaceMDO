function initialize_globals()
    % initialize_globals: Declare and initialize global variables

    global G M R_earth R_mars earth_orbital_data mars_result_data;

    % Define global constants
    G = 6.67430e-11;      % Gravitational constant (m^3/kg/s^2)
    M = 5.972e24;         % Mass of Earth (kg)
    R_earth = 6371e3;     % Radius of Earth (m)
    R_mars = 3389.5e3;    % Radius of Mars (m)
  
    % Load data from CSV files
    try
        earth_orbital_data = readmatrix('earth_result_table.csv');
        mars_result_data = readmatrix('mars_result_table.csv'); 

        % Display a confirmation message
        disp('Global variables and CSV data initialized successfully.');
    catch ME
        error('Error loading CSV files: %s', ME.message);
    end

end
