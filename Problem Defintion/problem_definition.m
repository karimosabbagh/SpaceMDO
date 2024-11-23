function PB = problem_definition

% The constant LAMBDA is also stored in PB
% This is a user defined parameter that is only used in the subsystem
% analysis and is independent from the MDO method.
% You can define as many such variables as you want.
PB.UserData.LAMBDA = 0;


lb = 0;
ub = 10;

%            Name     SP      CV      links    dim    lb     ub
%----------------------------------------------------------
PB.var{1}  = {'u_1'   1    false      'u_2'     1     lb     ub}; % <----- example of a variable shared between subsystems 1 and 2
PB.var{2}  = {'v'     1    false        []      1     lb     ub}; % <----- example of a local variable to subsystem 1
PB.var{3}  = {'a_1'   1     true      'a_2'     1     lb     ub}; % <----- example of an outgoing coupling variable from subsystem 1
PB.var{4}  = {'b_1'   1    false      'b_2'     1     lb     ub};
%----------------------------------------------------------
PB.var{5}   = {'u_2'  2    false      'u_1'     1     lb     ub}; % <----- example of a variable shared between subsystems 1 and 2
PB.var{6}   = {'w'    2    false        []      1     lb     ub}; % <----- example of a local variable to subsystem 2
PB.var{7}   = {'a_2'  2    false      'a_1'     1     lb     ub};
PB.var{8}   = {'b_2'  2     true      'b_1'     1     lb     ub}; % <----- example of an outgoing coupling variable from subsystem 2

% The objective function of sub-system index_main is considered as the 
% general objective function
PB.index_main = 1;
% Function to call to perform the subsystem analysis:
PB.analysis_file = 'example_subsystem_analysis';
PB.end_of_iter_file = 'example_display';