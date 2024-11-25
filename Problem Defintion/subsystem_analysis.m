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

switch subproblem_index
    case 1
        % Subproblem 1 - ORBITAL ESCAPE
        [m_SC, r_p1, V_SC_arrival, departure_date, arrival_date] = get_variable(x_DV,PB,'m_SC_e', 'r_p1', 'V_SC_arrival_e', 'departure_date_e', 'arrival_date_e');
        [delta_v_escape, V_SC_departure, S1_constraints] = orbital_escape_delta_v(m_SC, r_p1, V_SC_arrival, departure_date, arrival_date);
        obj = delta_v_escape;
        c_ineq(1) = S1_constraints(1);
        c_ineq(2) = S1_constraints(2);
        c_ineq(3) = S1_constraints(3);
        c_ineq(4) = S1_constraints(4);
        y = [delta_v_escape, V_SC_departure];
    case 2
        % Subproblem 2 - ORBITAL CAPTURE
        [u,w,a] = get_variable(x_DV,PB,'u_2','w','a_2');            
        b = 1/(u+LAMBDA)+1/(w+LAMBDA)+1/(a+LAMBDA);
        obj = 0;
        y = b;
        c_ineq = w+b-10;
    case 3
        % Subproblem 3 - SPACECRAFT & PROPELLANT MASS
        [u,w,a] = get_variable(x_DV,PB,'u_3','w','a_3');            
        b = 1/(u+LAMBDA)+1/(w+LAMBDA)+1/(a+LAMBDA);
        obj = 0;
        y = b;
        c_ineq = w+b-10;
     case 4
        % Subproblem 4 - PLANET COVERAGE
        [u,w,a] = get_variable(x_DV,PB,'u_4','w','a_4');            
        b = 1/(u+LAMBDA)+1/(w+LAMBDA)+1/(a+LAMBDA);
        obj = 0;
        y = b;
        c_ineq = w+b-10;
    otherwise
        error('unrecognized subproblem index')
end
