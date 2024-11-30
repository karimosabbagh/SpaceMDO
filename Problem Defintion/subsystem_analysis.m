function [obj,y,c_ineq] = subsystem_analysis(subproblem_index,x_DV,PB)

% Add subsystem paths
currentFilePath = fileparts(mfilename('fullpath'));
subsytems = fullfile(currentFilePath, '..', 'Subsystems');
addpath(subsytems);

% User defined constant
LAMBDA = PB.UserData.LAMBDA; 

% default output values
obj = 0;
y = [];
c_ineq = [];
%disp('subproblem index')
%disp(subproblem_index)
switch subproblem_index
    case 1
        % Subproblem 1 - ORBITAL ESCAPE
        [r_p1, V_SC_arrival, V_SC_departure, departure_date, arrival_date] = get_variable(x_DV,PB, 'r_p_e', 'V_SC_arrival_e', 'V_SC_departure_e', 'departure_date_e', 'arrival_date_e');
        [delta_v_escape, V_SC_departure, S1_constraints] = S1_orbital_escape(r_p1, V_SC_arrival, departure_date, arrival_date);
        obj = delta_v_escape;
        y = V_SC_departure;
        c_ineq = S1_constraints;

    case 2
        % Subproblem 2 - ORBITAL CAPTURE
        [r_p2, V_SC_departure, departure_date, arrival_date, e,delta_m_d] = get_variable(x_DV,PB, 'r_p_c', 'V_SC_departure_c', 'departure_date_c', 'arrival_date_c', 'e_c','m_prop_c');            
        [delta_v_capture, V_SC_arrival, S2_constraints] = S2_orbital_capture(e, delta_m_d, r_p2, V_SC_departure, departure_date, arrival_date);
        obj = delta_v_capture;
        c_ineq = S2_constraints;
        y = V_SC_arrival;

    case 3
        % Subproblem 3 - SPACECRAFT & PROPELLANT MASS
        [delta_v_escape,delta_v_arrival, departure_date, arrival_date, Isp, m_structure] = get_variable(x_DV,PB,'delta_v_escape_s','delta_v_capture_s','departure_date_s','arrival_date_s','Isp','m_structure_s');
        [m_prop_s,cost,S3_constraints] = S3_prop_struc_mass(delta_v_escape,delta_v_arrival, departure_date, arrival_date, Isp, m_structure);
        obj = cost;
        c_ineq = S3_constraints;
        y = m_prop_s;

     case 4
        % Subproblem 4 - PLANET COVERAGE
        [r_p, e, T_orbit, eta_center, eta_FOV_tilde, IFOV] = get_variable(x_DV, PB, 'r_p_p', 'e_p', 'T_orbit', 'eta_center', 'eta_FOV', 'IFOV');
        [percent_coverage, T_orbit, S4_constraints] = S4_planet_coverage(r_p, e, eta_center, eta_FOV_tilde, IFOV);
        obj = percent_coverage;
        % y = T_orbit;
        y = [];
        c_ineq = S4_constraints;
    otherwise
        error('unrecognized subproblem index')
end
