%-------------------------------------------------------------------------------------%
%  NoHiMDO                                                                            %
%                                                                                     %
%  A solver for Multi-Disciplinary Optimization, based on Non-Hierarchical Analytical %
%  Target Cascading                                                                   %
%  Version 3.0.0                                                                      %
%                                                                                     %
%  Copyright (C) 2012-2019  Bastien Talgorn - McGill University, Montreal             %
%                                                                                     %
%  Author: Bastien Talgorn                                                            %
%  email: bastientalgorn@fastmail.com                                                 %
%                                                                                     %
%  This program is free software: you can redistribute it and/or modify it under the  %
%  terms of the GNU Lesser General Public License as published by the Free Software   %
%  Foundation, either version 3 of the License, or (at your option) any later         %
%  version.                                                                           %
%                                                                                     %
%  This program is distributed in the hope that it will be useful, but WITHOUT ANY    %
%  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A    %
%  PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.   %
%                                                                                     %
%  You should have received a copy of the GNU Lesser General Public License along     %
%  with this program. If not, see <http://www.gnu.org/licenses/>.                     %
%                                                                                     %
%  You can find information on NoHiMDO at https://github.com/bastientalgorn/NoHiMDO   %
%-------------------------------------------------------------------------------------%

clc
close all
clear all
addpath NoHIMDO_solver
addpath Setup
addpath Problem_Definition
addpath Analysis

% Create a timestamped folder path under 'Experiments'
timestamp = datestr(now, 'yyyy-mm-dd_HH-MM-SS');  % Correct timestamp format
folder_path = fullfile('Experiments', timestamp);  % Create path under 'Experiments'

% Create the folder if it doesn't exist
if ~exist(folder_path, 'dir')
    mkdir(folder_path);
end

% Initialize variables
%global G M R_earth R_mars earth_orbital_data mars_result_data;
%x = [];
initialize_globals();
[x, initial_guess] = generate_initial_guess();

% Save the initial guess structure to the specified folder 
IG_filename = 'initial_guess.mat'; % Add file extension to save as .mat
save(fullfile(folder_path, IG_filename), 'initial_guess');


% NiHiMDO Parameters
NoHi_options.display = true;
NoHi_options.w_scheme = 'median';
NoHi_options.constraints_cv = true;
NoHi_options.realistic_obj = false;
NoHi_options.noprogress_stop = 100;
NoHi_options.NI = 1000;
NoHi_options.NO = 100;
NoHi_options.beta = 1.3;
NoHi_options.gamma = 0.5;
NoHi_options.w0 = 1;
NoHi_options.x0 = x;
NoHi_options.inc_stop = 1e-12;
NoHi_options.tol = 1e-6;
NoHi_options.nb_proc = 1;
NoHi_options.save_subproblems = false;
NoHi_options.solver = 'interior-point'; % options are : 'mads','sqp','interior-point','active-set','trust-region-reflective'
NoHi_options.solver_display = true;

% Conditional parallel pool initialization
if NoHi_options.nb_proc > 1
    % Start parallel pool if not already started
    if isempty(gcp('nocreate'))
        parpool(NoHi_options.nb_proc); % Start a pool with the specified number of workers
    end
    % Broadcast global variables to workers
    parfevalOnAll(@initialize_globals, 0);
else
    % Initialize global variables for single-threaded execution
    initialize_globals();
end


PB = problem_definition;
output = NoHiSolver(PB,NoHi_options);

% save output.mat output PB NoHi_options
% Save to appropriate timestamp folder
save(fullfile(folder_path, 'output.mat'), 'output', 'PB', 'NoHi_options');

% save results in neat formatted manner
results_table = format_results(output, PB);
save(fullfile(folder_path, 'results.mat'), 'results_table');


% Evaluate Subsystems
[results, evaluations] =  evaluate_subsystems(output.x, PB);
save(fullfile(folder_path, 'evaluations.mat'), 'evaluations');