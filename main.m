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
addpath Setup\
addpath 'Problem Defintion'\

% Initialize variables
%global G M R_earth R_mars earth_orbital_data mars_result_data;
initialize_globals();

% NiHiMDO Parameters
NoHi_options.display = true;
NoHi_options.w_scheme = 'median';
NoHi_options.constraints_cv = true;
NoHi_options.realistic_obj = false;
NoHi_options.noprogress_stop = 100;
NoHi_options.NI = 1000;
NoHi_options.NO = 200;
NoHi_options.beta = 1.3;
NoHi_options.gamma = 0.5;
NoHi_options.w0 = 1;
NoHi_options.x0 = [];
NoHi_options.inc_stop = 1e-12;
NoHi_options.tol = 1e-12;
NoHi_options.nb_proc = 1;
NoHi_options.save_subproblems = false;
NoHi_options.solver = 'mads'; % options are : 'mads','sqp','interior-point','active-set','trust-region-reflective'
NoHi_options.solver_display = false;


PB = problem_definition;
output = NoHiSolver(PB,NoHi_options);
output;

save output.mat output PB NoHi_options
