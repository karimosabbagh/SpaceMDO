function PB = problem_definition

% The constant LAMBDA is also stored in PB
% This is a user defined parameter that is only used in the subsystem
% analysis and is independent from the MDO method.
% You can define as many such variables as you want.
PB.UserData.LAMBDA = 0;

global R_mars R_earth start_date end_date;

lb = 0;
ub = inf;


%              Name                   SP    CV    links       dim   lb               ub
%
% ORBITAL ESCAPE
%----------------------------------------------------------
PB.var{1}   = {'delta_v_escape_e'     1    true     []         1    lb                80e3};
% PB.var{2}   = {'m_SC_e'             1    false    []         1    lb                ub};  
PB.var{2}   = {'r_p_e'                1    false    []         1    R_earth + 160e3   R_earth + 2000e3}; 
PB.var{3}   = {'departure_date_e'     1    false    []         1    start_date        end_date};
PB.var{4}   = {'arrival_date_e'       1    false    []         1    start_date        end_date};
PB.var{5}   = {'V_SC_departure_e'     1    true     [13]         1    lb                80e3};
PB.var{6}   = {'V_SC_arrival_e'       1    false    [14]         1    lb                80e3};

% ORBITAL CAPTURE
%----------------------------------------------------------
PB.var{7}   = {'delta_v_capture_c'    2    true     []         1    lb                80e3};
% PB.var{9}   = {'m_SC_c'             2    false    []         1    lb                80e3}; 
PB.var{8}  =  {'m_prop_c'             2    false    []         1    lb                80e3}; 
PB.var{9}  =  {'r_p_c'                2    false    'r_p_p'         1    R_mars + 100e3    170 * R_mars}; 
PB.var{10}  = {'e_c'                  2    false    []         1    lb                1};
PB.var{11}  = {'departure_date_c'     2    false    []         1    start_date        end_date};
PB.var{12}  = {'arrival_date_c'       2    false    []         1    start_date        end_date};
PB.var{13}  = {'V_SC_departure_c'     2    false    [5]         1    lb                80e3};
PB.var{14}  = {'V_SC_arrival_c'       2    true     [6]         1    lb                80e3};

% SPACECRAFT & PROPELLANT MASS
%----------------------------------------------------------
PB.var{15}   = {'cost'                3    false    []        1    lb                 80e3}; 
% PB.var{18}   = {'m_SC_s'            3    true     []        1    lb                 3000};
%PB.var{16}   = {'m_prop_s'            3    false     []       1    lb                 1000};
PB.var{17}   = {'m_structure_s'       3    false    []        1    lb                 1000};
PB.var{18}   = {'Isp'                 3    false    []        1    1                  600}; 
PB.var{19}   = {'delta_v_escape_s'    3    false    []        1    lb                 80e3};
PB.var{20}   = {'delta_v_capture_s'   3    false    []        1    lb                 80e3};
PB.var{21}   = {'departure_date_s'    3    false    []        1    start_date         end_date};
PB.var{22}   = {'arrival_date_s'      3    false    []        1    start_date         end_date};
PB.var{23}   = {'m_payload'           3    false    []        1    lb                 1000};

% PLANET COVERAGE
%----------------------------------------------------------
PB.var{24}   = {'r_p_p'               4    false    'r_p_e'        1    R_mars             170 * R_mars  };
PB.var{25}   = {'e_p'                 4    false    []        1    lb                 1};
PB.var{26}   = {'T_orbit'             4    false    []        1    1e3                10e3};  
PB.var{27}   = {'eta_center'          4    false    []        1    lb                 deg2rad(45)};
PB.var{28}   = {'eta_FOV'             4    false    []        1    lb                 deg2rad(35)};
PB.var{29}   = {'IFOV'                4    false    []        1    1e-3               10};



% %              Name                   SP    CV    links       dim   lb               ub
% %
% % ORBITAL ESCAPE
% %----------------------------------------------------------
% PB.var{1}   = {'delta_v_escape_e'     1    true     [19]       1    lb                80e3};
% % PB.var{2}   = {'m_SC_e'             1    false    [9,18]     1    lb                ub};  
% PB.var{2}   = {'r_p_e'                1    false    []         1    R_earth + 160e3   R_earth + 2000e3}; 
% PB.var{3}   = {'departure_date_e'     1    false    [11,21]    1    start_date        end_date};
% PB.var{4}   = {'arrival_date_e'       1    false    [12,22]    1    start_date        end_date};
% PB.var{5}   = {'V_SC_departure_e'     1    true     [13]       1    lb                80e3};
% PB.var{6}   = {'V_SC_arrival_e'       1    false    [14]       1    lb                80e3};
% 
% % ORBITAL CAPTURE
% %----------------------------------------------------------
% PB.var{7}   = {'delta_v_capture_c'    2    true     [20]       1    lb                80e3};
% % PB.var{9}   = {'m_SC_c'             2    false    [2,18]     1    lb                80e3}; 
% PB.var{8}  =  {'m_prop_c'             2    false    [16]       1    lb                80e3}; 
% PB.var{9}  =  {'r_p_c'                2    false    [24]       1    R_mars + 100e3    170 * R_mars}; 
% PB.var{10}  = {'e_c'                  2    false    [25]         1    lb                1};
% PB.var{11}  = {'departure_date_c'     2    false    [3,21]     1    start_date        end_date};
% PB.var{12}  = {'arrival_date_c'       2    false    [4,22]     1    start_date        end_date};
% PB.var{13}  = {'V_SC_departure_c'     2    false    [5]        1    lb                80e3};
% PB.var{14}  = {'V_SC_arrival_c'       2    true     [6]        1    lb                80e3};
% 
% % SPACECRAFT & PROPELLANT MASS
% %----------------------------------------------------------
% PB.var{15}   = {'cost'                3    false    []        1    lb                 80e3}; 
% % PB.var{18}   = {'m_SC_s'            3    true     [2,9]     1    lb                 3000};
% PB.var{16}   = {'m_prop_s'            3    true     [8]      1    lb                 1000};
% PB.var{17}   = {'m_structure_s'       3    false    []        1    lb                 1000};
% PB.var{18}   = {'Isp'                 3    false    []        1    1                  600}; 
% PB.var{19}   = {'delta_v_escape_s'    3    false    [1]       1    lb                 80e3};
% PB.var{20}   = {'delta_v_capture_s'   3    false    [7]       1    lb                 80e3};
% PB.var{21}   = {'departure_date_s'    3    false    [3,11]    1    start_date         end_date};
% PB.var{22}   = {'arrival_date_s'      3    false    [4,12]    1    start_date         end_date};
% PB.var{23}   = {'m_payload'           3    false    []        1    lb                 1000};
% 
% 
% % PLANET COVERAGE
% %----------------------------------------------------------
% PB.var{24}   = {'r_p_p'               4    false    'r_p_c'   1    R_mars             170 * R_mars  };
% PB.var{25}   = {'e_p'                 4    false    'e_c'     1    lb                 1};
% PB.var{26}   = {'T_orbit'             4    false    []        1    1e3                10e3};  
% PB.var{27}   = {'eta_center'          4    false    []        1    lb                 deg2rad(45)};
% PB.var{28}   = {'eta_FOV'             4    false    []        1    lb                 deg2rad(35)};
% PB.var{29}   = {'IFOV'                4    false    []        1    1e-3               10};
% 

% The objective function of sub-system index_main is considered as the 
% general objective function
PB.index_main = 3;
% Function to call to perform the subsystem analysis:
PB.analysis_file = 'subsystem_analysis';
PB.end_of_iter_file = 'display';


% PLANET COVERAGE
%----------------------------------------------------------
%PB.var{29}   = {'u_2'  4    false      []     1     lb     ub}; % <----- example of a variable shared between subsystems 1 and 2
%PB.var{30}   = {'w'    4    false        []      1     lb     ub}; % <----- example of a local variable to subsystem 2
%PB.var{31}   = {'a_4'  4    false      'a_1'     1     lb     ub};
%PB.var{32}   = {'b_4'  4     true      'b_1'     1     lb     ub}; % <----- example of an outgoing coupling variable from subsystem 2