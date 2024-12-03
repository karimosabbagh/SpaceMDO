function PB = problem_definition

% The constant LAMBDA is also stored in PB
% This is a user defined parameter that is only used in the subsystem
% analysis and is independent from the MDO method.
% You can define as many such variables as you want.
PB.UserData.LAMBDA = 0;

global R_mars R_earth;

lb = 0;
ub = inf;

%              Name                   SP    CV     links                     dim  lb               ub
%
% ORBITAL ESCAPE
%----------------------------------------------------------
PB.var{1}   = {'delta_v_escape_e'     1    true     'delta_v_escape_s'       1    lb                80e3};
PB.var{2}   = {'r_p_e'                1    false    []                       1    (R_earth+160e3)   (R_earth+2000e3)}; 
%PB.var{3}   = {'departure_date_e'     1    false    [11,20]    1    start_date        (end_date - 100)};
%PB.var{4}   = {'arrival_date_e'       1    false    [12,21]    1    (start_date + 100)        end_date};
PB.var{3}   = {'V_SC_departure_e'     1    true     'V_SC_departure_c'       1    lb                80e3};
PB.var{4}   = {'V_SC_arrival_e'       1    false    'V_SC_arrival_c'         1    lb                80e3};

% ORBITAL CAPTURE
%----------------------------------------------------------
PB.var{5}  =  {'delta_v_capture_c'    2    true     'delta_v_capture_s'      1    lb                80e3};
PB.var{6}  =  {'m_prop_c'             2    false    'm_prop_s'               1    lb                3000}; 
PB.var{7}  =  {'r_p_c'                2    false    'r_p_p'                  1    (R_mars+100e3)    (1000e3+R_mars)}; 
PB.var{8}  =  {'e_c'                  2    false    'e_p'                    1    lb                0.9};
%PB.var{11}  = {'departure_date_c'     2    false    [3,20]     1    start_date        (end_date - 100)};
%PB.var{12}  = {'arrival_date_c'       2    false    [4,21]     1    (start_date + 100)        end_date};
PB.var{9}  =  {'V_SC_departure_c'     2    false    'V_SC_departure_e'       1    lb                80e3};
PB.var{10} =  {'V_SC_arrival_c'       2    true     'V_SC_arrival_e'         1    lb                80e3};

% SPACECRAFT & PROPELLANT MASS
%----------------------------------------------------------
PB.var{11}   = {'m_prop_s'            3    true     'm_prop_c'               1    lb                3000}; 
PB.var{12}   = {'m_structure_s'       3    false    []                       1    lb                1500};
PB.var{13}   = {'Isp'                 3    false    []                       1    150               600}; 
PB.var{14}   = {'delta_v_escape_s'    3    false    'delta_v_escape_e'       1    lb                80e3};
PB.var{15}   = {'delta_v_capture_s'   3    false    'delta_v_capture_c'      1    lb                80e3};
%PB.var{20}   = {'departure_date_s'    3    false    [3,11]    1    start_date         (end_date - 100)};
%PB.var{21}   = {'arrival_date_s'      3    false    [4,12]    1    (start_date + 100)         end_date};
PB.var{16}   = {'eta_FOV_s'           3    false    'eta_FOV_p'              1    lb                deg2rad(35)};
PB.var{17}   = {'IFOV_s'              3    false    'IFOV_p'                 1    1e-3              10};


% PLANET COVERAGE
%----------------------------------------------------------
PB.var{18}   = {'r_p_p'               4    false    'r_p_c'                  1    (R_mars+100e3)   (1000e3+R_mars)};
PB.var{19}   = {'e_p'                 4    false    'e_c'                    1    lb                0.9};
%PB.var{20}   = {'T_orbit'             4    false    []                       1    1e3               10e4};  
PB.var{20}   = {'eta_center'          4    false    []                       1    lb                deg2rad(45)};
PB.var{21}   = {'eta_FOV_p'           4    false    'eta_FOV_s'              1    lb                deg2rad(35)};
PB.var{22}   = {'IFOV_p'              4    false    'IFOV_s'                 1    1e-3              10};


% The objective function of sub-system index_main is considered as the 
% general objective function
PB.index_main = 3;
% Function to call to perform the subsystem analysis:
PB.analysis_file = 'subsystem_analysis';
PB.end_of_iter_file = 'display_convergence';