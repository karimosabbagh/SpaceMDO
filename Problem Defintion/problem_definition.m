function PB = problem_definition

% The constant LAMBDA is also stored in PB
% This is a user defined parameter that is only used in the subsystem
% analysis and is independent from the MDO method.
% You can define as many such variables as you want.
PB.UserData.LAMBDA = 0;

global date_range R_mars;

%lb = 0;
%ub = 10;



%              Name                   SP    CV    links  dim   lb    ub
%
% ORBITAL ESCAPE
%----------------------------------------------------------
PB.var{1}   = {'delta_v_escape'       1     true    []      1    0    ub};
PB.var{2}   = {'m_SC_e'               1    false    [10]    1    lb    ub};  % from Spacecraft & Propellant Mass
PB.var{3}   = {'r_p1'                 1    false    []      1    160000    2000000}; 
PB.var{4}   = {'departure_date_e'     1    false    [14]    1    date_range(1)    date_range(end)};
PB.var{5}   = {'arrival_date_e'       1    false    [15]    1    date_range(1)    date_range(end)};
PB.var{6}   = {'tof_e'                1    false    [16]    1    128    500};
PB.var{7}   = {'V_SC_departure_e'     1     true    [17]    1    0    ub};
PB.var{8}   = {'V_SC_arrival_e'       1    false    [18]    1    0    ub};

% ORBITAL CAPTURE
%----------------------------------------------------------
PB.var{9}   = {'delta_v_capture'      2     true    []     1    0    ub};
PB.var{10}  = {'m_SC_c'               2    false    [2]    1    lb    ub};  % from Spacecraft & Propellant Mass
PB.var{11}  = {'m_prop_c'             2    false    []     1    lb    ub};  % from Spacecraft & Propellant Mass
PB.var{12}  = {'r_p2'                 2    false    []     1    100000    320000}; 
PB.var{13}  = {'e_2'                  2    false    []     1    0    1};
PB.var{14}  = {'departure_date_c'     2    false    [4]    1    date_range(1)    date_range(end)};
PB.var{15}  = {'arrival_date_c'       2    false    [5]    1    date_range(1)    date_range(end)};
PB.var{16}  = {'tof_c'                2    false    [6]    1    128    500};
PB.var{17}  = {'V_SC_departure_c'     2    false    [7]    1    0    ub};
PB.var{18}  = {'V_SC_arrival_c'       2     true    [8]    1    0    ub};

% SPACECRAFT & PROPELLANT MASS
%----------------------------------------------------------
PB.var{25}   = {'u_3'                 3    false      'u_1'     1     lb     ub}; % <----- example of a variable shared between subsystems 1 and 2
PB.var{26}   = {'w'                   3    false        []      1     lb     ub}; % <----- example of a local variable to subsystem 2
PB.var{27}   = {'a_3'                 3    false      'a_1'     1     lb     ub};
PB.var{28}   = {'b_3'                 3     true      'b_1'     1     lb     ub}; % <----- example of an outgoing coupling variable from subsystem 2

% PLANET COVERAGE
%----------------------------------------------------------
PB.var{29}   = {'r_p3'                4    false     [13]     1     R_mars   170 * R_mars  }; 
PB.var{31}   = {'T_orbit_4'           4    true     'T_orbit_3'     1e3     inf};  % Outgoing coupling variable to S3 
PB.var{32}   = {'eta_center'          4    false     []       1     0     deg2rad(45)};
PB.var{33}   = {'eta_FOV'             4    false     []       1     0     deg2rad(35)};
PB.var{34}   = {'IFOV'                4    false     []       1     1e-3  10};

% The objective function of sub-system index_main is considered as the 
% general objective function
PB.index_main = 1;
% Function to call to perform the subsystem analysis:
PB.analysis_file = 'subsystem_analysis';
PB.end_of_iter_file = 'display';