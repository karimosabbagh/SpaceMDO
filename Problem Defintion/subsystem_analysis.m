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
        % Subproblem 1
        [u,v,b] = get_variable(x_DV,PB,'u_1','v','b_1');  
        a = log(u)+log(b)+log(v);
        obj = u+v+a+b; 
        y = a;
        c_ineq = [];
    case 2
        % Subproblem 2
        [u,w,a] = get_variable(x_DV,PB,'u_2','w','a_2');            
        b = 1/(u+LAMBDA)+1/(w+LAMBDA)+1/(a+LAMBDA);
        obj = 0;
        y = b;
        c_ineq = w+b-10;
    otherwise
        error('unrecognized subproblem index')
end
