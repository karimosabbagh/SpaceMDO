function [results, evaluations] =  evaluate_subsystems(x, PB)
        
    % Add subsystem paths
    currentFilePath = fileparts(mfilename('fullpath'));
    subsytems = fullfile(currentFilePath, '..', 'Subsystems');
    addpath(subsytems);

    % Get variable names
    numVars = numel(PB.var);
    results = struct(); 

    for i = 1:numVars
        varName = PB.var{i}{1}; % Get the variable name
        results.(varName) = x(i); % Dynamically add fields to the struct
    end
    disp(results);
    evaluations = struct();
    
    % Subsystem 1
    [evaluations.delta_v_escape, evaluations.V_SC_departure, evaluations.S1_constraints] ... 
        = S1_orbital_escape(results.r_p_e, results.V_SC_arrival_e, results.departure_date_e, results.arrival_date_e);
    
    % Subsystem 2
    [evaluations.delta_v_capture, evaluations.V_SC_arrival, evaluations.S2_constraints]   ...
        = S2_orbital_capture(results.e_c, results.m_prop_c, results.r_p_c, results.V_SC_departure_c, results.departure_date_c, results.arrival_date_c);
    
    % Subsystem 3
    [evaluations.m_prop_s, evaluations.cost, evaluations.S3_constraints] = ...
        S3_prop_struc_mass(results.delta_v_escape_s, results.delta_v_capture_s, results.departure_date_s, results.arrival_date_s, results.Isp, results.m_structure_s);

    % Subsystem 4
    [evaluations.percent_coverage, evaluations.T_orbit, evaluations.S4_constraints] = ...
        S4_planet_coverage(results.r_p_p, results.e_p, results.eta_center, results.eta_FOV, results.IFOV);

    % Display the results struct for verification
    disp(evaluations);
end