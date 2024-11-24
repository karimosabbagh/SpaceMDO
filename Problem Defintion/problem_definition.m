function PB = problem_definition

% The constant LAMBDA is also stored in PB
% This is a user defined parameter that is only used in the subsystem
% analysis and is independent from the MDO method.
% You can define as many such variables as you want.
PB.UserData.LAMBDA = 0;

global date_range;

%lb = 0;
%ub = 10;



%              Name                   SP    CV    links  dim   lb    ub
%
% ORBITAL ESCAPE
%----------------------------------------------------------
PB.var{1}   = {'delta_v_escape'       1     true    []      1    0    ub};
PB.var{2}   = {'m_SC_e'               1    false    [13]    1    lb    ub};  % from Spacecraft & Propellant Mass
PB.var{3}   = {'r_p1'                 1    false    []      1    160000    2000000}; 
PB.var{4}   = {'departure_date_e'     1    false    [17]    1    date_range(1)    date_range(end)};
PB.var{5}   = {'arrival_date_e'       1    false    [18]    1    date_range(1)    date_range(end)};
PB.var{6}   = {'tof_e'                1    false    [19]    1    128    500};
PB.var{7}   = {'V_SC_departure_e'     1     true    [20]    1    0    ub};
PB.var{8}   = {'V_SC_arrival_e'       1    false    [21]    1    0    ub};
%PB.var{9}   = {'V_Earth_departure'    1    false    []    1    29291.3527030927    30286.3195731899};         % not sure if this should count, as it depends on departure date. They are linked.
%PB.var{10}  = {'R_Earth_departure'    1    false    []    1    147099571771.6730042    152096240315.759979};  % not sure if this should count, as it depends on departure date. They are linked.
%PB.var{11}  = {'R_Mars_arrival'       1    false    []    1    206635022887.496    249232882114.895};         % not sure if this should count, as it depends on arrival date. They are linked.

% ORBITAL CAPTURE
%----------------------------------------------------------
PB.var{12}  = {'delta_v_capture'      2     true    []     1    0    ub};
PB.var{13}  = {'m_SC_c'               2    false    [2]    1    lb    ub};  % from Spacecraft & Propellant Mass
PB.var{14}  = {'m_prop_c'             2    false    []     1    lb    ub};  % from Spacecraft & Propellant Mass
PB.var{15}  = {'r_p2'                 2    false    []     1    100000    320000}; 
PB.var{16}  = {'e'                    2    false    []     1    100000    320000};
PB.var{17}  = {'departure_date_c'     2    false    [4]    1    date_range(1)    date_range(end)};
PB.var{18}  = {'arrival_date_c'       2    false    [5]    1    date_range(1)    date_range(end)};
PB.var{19}  = {'tof_c'                2    false    [6]    1    128    500};
PB.var{20}  = {'V_SC_departure_c'     2    false    [7]    1    0    ub};
PB.var{21}  = {'V_SC_arrival_c'       2     true    [8]    1    0    ub};
%PB.var{22}  = {'V_Mars_departure'     2    false    []    1    21971.0215481833    26500.3632618652};         % not sure if this should count, as it depends on arrival date. They are linked.
%PB.var{23}  = {'R_Mars_arrival'       2    false    []    1    206635022887.496    249232882114.895};         % not sure if this should count, as it depends on arrival date. They are linked.
%PB.var{24}  = {'R_Earth_departure'    2    false    []    1    147099571771.6730042    152096240315.759979};  % not sure if this should count, as it depends on departure date. They are linked.

% SPACECRAFT & PROPELLANT MASS
%----------------------------------------------------------
PB.var{25}   = {'u_3'  3    false      'u_1'     1     lb     ub}; % <----- example of a variable shared between subsystems 1 and 2
PB.var{26}   = {'w'    3    false        []      1     lb     ub}; % <----- example of a local variable to subsystem 2
PB.var{27}   = {'a_3'  3    false      'a_1'     1     lb     ub};
PB.var{28}   = {'b_3'  3     true      'b_1'     1     lb     ub}; % <----- example of an outgoing coupling variable from subsystem 2

% PLANET COVERAGE
%----------------------------------------------------------
PB.var{29}   = {'u_2'  4    false      'u_1'     1     lb     ub}; % <----- example of a variable shared between subsystems 1 and 2
PB.var{30}   = {'w'    4    false        []      1     lb     ub}; % <----- example of a local variable to subsystem 2
PB.var{31}   = {'a_4'  4    false      'a_1'     1     lb     ub};
PB.var{32}   = {'b_4'  4     true      'b_1'     1     lb     ub}; % <----- example of an outgoing coupling variable from subsystem 2

% The objective function of sub-system index_main is considered as the 
% general objective function
PB.index_main = 1;
% Function to call to perform the subsystem analysis:
PB.analysis_file = 'example_subsystem_analysis';
PB.end_of_iter_file = 'example_display';