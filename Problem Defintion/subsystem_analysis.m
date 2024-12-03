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
        [r_p1, V_SC_arrival] = get_variable(x_DV,PB, 'r_p_e', 'V_SC_arrival_e');
        [delta_v_escape, V_SC_departure, V_infinity_D] = S1_orbital_escape(r_p1, V_SC_arrival);
        obj = delta_v_escape;
        y = [V_SC_departure, delta_v_escape];
        c_ineq = V_infinity_D;

    case 2
        % Subproblem 2 - ORBITAL CAPTURE
        [r_p2, V_SC_departure, e, m_prop] = get_variable(x_DV,PB, 'r_p_c', 'V_SC_departure_c', 'e_c','m_prop_c');            
        [delta_v_capture, V_SC_arrival, V_infinity_A] = S2_orbital_capture(e, m_prop, r_p2, V_SC_departure);
        obj = delta_v_capture;
        c_ineq = V_infinity_A;
        y = [V_SC_arrival, delta_v_capture];

    case 3
        % Subproblem 3 - SPACECRAFT & PROPELLANT MASS
        [delta_v_escape, delta_v_arrival, Isp, m_structure, eta_FOV, IFOV] = get_variable(x_DV,PB,'delta_v_escape_s','delta_v_capture_s','Isp','m_structure_s','eta_FOV_s','IFOV_s');
        [m_prop, cost, shield] = S3_prop_struc_mass(delta_v_escape, delta_v_arrival, Isp, m_structure, eta_FOV, IFOV);
        obj = cost;
        c_ineq = shield;
        y = m_prop;

     case 4
        % Subproblem 4 - PLANET COVERAGE
        [r_p, e, eta_center, eta_FOV_tilde, IFOV] = get_variable(x_DV, PB, 'r_p_p', 'e_p', 'eta_center', 'eta_FOV_p', 'IFOV_p');
        [percent_coverage, S4_constraints] = S4_planet_coverage(r_p, e, eta_center, eta_FOV_tilde, IFOV);
        obj = percent_coverage;
        y = [];
        c_ineq = S4_constraints;
    otherwise
        error('unrecognized subproblem index')
end
