function [results, evaluations] =  evaluate_subsystems(x, PB)
        
    % Add subsystem paths
    currentFilePath = fileparts(mfilename('fullpath'));
    subsytems = fullfile(currentFilePath, '..', 'Subsystems');
    addpath(subsytems);

    start_date = datetime('2024-11-21', 'InputFormat', 'yyyy-MM-dd');
    departure_date = datenum(start_date);
    departure_date_e = fix(departure_date);
    departure_date_c = departure_date_e;
    departure_date_s = departure_date_e;
    end_date = datetime('2025-05-20', 'InputFormat', 'yyyy-MM-dd');
    arrival_date = datenum(end_date);
    arrival_date_e = fix(arrival_date);
    arrival_date_c = arrival_date_e;
    arrival_date_s = arrival_date_e;

    % Get variable names
    numVars = numel(PB.var);
    results = struct(); 

    for i = 1:numVars
        varName = PB.var{i}{1}; % Get the variable name
        results.(varName) = x(i); % Dynamically add fields to the struct
    end

    evaluations = struct();
    
    % Subsystem 1
    [evaluations.delta_v_escape, evaluations.V_SC_departure, evaluations.S1_constraints] ... 
        = S1_orbital_escape(results.r_p_e, results.V_SC_arrival_e, departure_date_e, arrival_date_e);
    
    % Subsystem 2
    [evaluations.delta_v_capture, evaluations.V_SC_arrival, evaluations.S2_constraints]   ...
        = S2_orbital_capture(results.e_c, results.r_p_c, results.V_SC_departure_c, departure_date_c, arrival_date_c);
    
    % Subsystem 3
    [evaluations.m_prop_s, evaluations.cost, evaluations.S3_constraints] = ...
        S3_prop_struc_mass(results.delta_v_escape_s, results.delta_v_capture_s, departure_date_s, arrival_date_s, results.Isp, results.m_structure_s, results.eta_FOV_s, results.IFOV_s);

    % Subsystem 4
    [evaluations.percent_coverage, evaluations.S4_constraints] = ...
        S4_planet_coverage(results.r_p_p, results.e_p, results.eta_center, results.eta_FOV_p, results.IFOV_p);

    % Display the results struct for verification
    disp(evaluations);
end