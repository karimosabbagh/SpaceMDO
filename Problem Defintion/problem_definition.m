function PB = problem_definition

% The constant LAMBDA is also stored in PB
% This is a user defined parameter that is only used in the subsystem
% analysis and is independent from the MDO method.
% You can define as many such variables as you want.
PB.UserData.LAMBDA = 0;

global date_range R_mars;

lb = 0;
ub = inf;



%              Name                   SP    CV    links       dim   lb               ub
%
% ORBITAL ESCAPE
%----------------------------------------------------------
PB.var{1}   = {'delta_v_escape'       1    true     [24]         1    lb                ub};
PB.var{2}   = {'m_SC_e'               1    false    [10,20]    1    lb                ub};  % from Spacecraft & Propellant Mass
PB.var{3}   = {'r_p_e'                1   false     []         1    R_earth + 160e3   R_earth + 2000e3}; 
PB.var{4}   = {'departure_date_e'     1    false    [14]       1    date_range(1)     date_range(end)};
PB.var{5}   = {'arrival_date_e'       1    false    [15]       1    date_range(1)     date_range(end)};
PB.var{6}   = {'tof_e'                1    true     [16, 26]   1    128               500};
PB.var{7}   = {'V_SC_departure_e'     1    true     [17]       1    lb                ub};
PB.var{8}   = {'V_SC_arrival_e'       1    false    [18]       1    lb                ub};

% ORBITAL CAPTURE
%----------------------------------------------------------
PB.var{9}   = {'delta_v_capture'      2    true     [25]        1    lb                 ub};
PB.var{10}  = {'m_SC_c'               2    false    [2,20]    1    lb                 ub};  % from Spacecraft & Propellant Mass
PB.var{11}  = {'m_prop_c'             2    false    [21]      1    lb                 ub};  % from Spacecraft & Propellant Mass
PB.var{12}  = {'r_p_c'                2    false    [27]      1    R_mars + 100e3     170 * R_mars}; 
PB.var{13}  = {'e_c'                  2    false    []        1    lb                 1};
PB.var{14}  = {'departure_date_c'     2    false    [4]       1    date_range(1)      date_range(end)};
PB.var{15}  = {'arrival_date_c'       2    false    [5]       1    date_range(1)      date_range(end)};
PB.var{16}  = {'tof_c'                2    false    [6, 26]       1    128                500};
PB.var{17}  = {'V_SC_departure_c'     2    false    [7]       1    lb                 ub};
PB.var{18}  = {'V_SC_arrival_c'       2    true     [8]       1    lb                 ub};

% SPACECRAFT & PROPELLANT MASS
%----------------------------------------------------------
PB.var{19}   = {'cost'                3    false    []        1    lb                 ub}; 
PB.var{20}   = {'m_SC_s'              3    true     [2,10]    1    lb                 3000}; 
PB.var{21}   = {'m_prop_s'            3    true     [11]      1    lb                 1000};
PB.var{22}   = {'m_structure_s'       3    false    []        1    lb                 1000};
PB.var{23}   = {'Isp'                 3    false    []        1    1                  600}; 
PB.var{24}   = {'delta_v_escape_e'    3    false    [1]       1    lb                 ub};
PB.var{25}   = {'delta_v_arrival_e'   3    false    [9]       1    lb                 ub};
PB.var{26}   = {'tof_s'               3    false    [6,16]    1    128                500};

% PLANET COVERAGE
%----------------------------------------------------------
PB.var{27}   = {'r_p_p'               4    false    [12]      1    R_mars             170 * R_mars  };
PB.var{28}   = {'e_p'                 2    false    [13]      1    lb                 1};
PB.var{29}   = {'T_orbit'             4    false    []        1    1e3                ub};  
PB.var{30}   = {'eta_center'          4    false    []        1    lb                 deg2rad(45)};
PB.var{31}   = {'eta_FOV'             4    false    []        1    lb                 deg2rad(35)};
PB.var{32}   = {'IFOV'                4    false    []        1    1e-3               10};


% The objective function of sub-system index_main is considered as the 
% general objective function
PB.index_main = 1;
% Function to call to perform the subsystem analysis:
PB.analysis_file = 'subsystem_analysis';
PB.end_of_iter_file = 'display';


% PLANET COVERAGE
%----------------------------------------------------------
%PB.var{29}   = {'u_2'  4    false      []     1     lb     ub}; % <----- example of a variable shared between subsystems 1 and 2
%PB.var{30}   = {'w'    4    false        []      1     lb     ub}; % <----- example of a local variable to subsystem 2
%PB.var{31}   = {'a_4'  4    false      'a_1'     1     lb     ub};
%PB.var{32}   = {'b_4'  4     true      'b_1'     1     lb     ub}; % <----- example of an outgoing coupling variable from subsystem 2