function [sys, data] = SimEngine2D(model_file, analysis_file)
% SIMENGINE2D - dynamic and kinematic analysis of planar mechanisms.
%   [SYS, DATA] = SimEngine2D(MODEL_FILENAME, ANALYSIS_FILENAME) performs
%   the analysis specified in the ACF file ANALYSIS_FILENAME using the
%   model described in the ADM file MODEL_FILENAME.
%
%   On completion, SYS is an object of type MBsys created from the
%   information in the ADM file (see MBsys for deatils), while DATA is a 
%   structure with simulation results (see KINEMATICS and DYNAMICS for 
%   details).
%
%   See also MBSYS, KINEMATICS, DYNAMICS

%% Create the MBsystem using the definition in the ADM file.
tic;
sys = MBsys(model_file);
fprintf('Time to read model: %f\n', toc);

%% Parse the ACF file to extract analysis information.
an = loadjson(analysis_file);

tic;
if strcmpi(an.simulation, 'kinematics')
    data = kinematics(sys, 0, an.tend, an.stepSize, an.outputSteps);
elseif strcmpi(an.simulation, 'dynamics')
    data = dynamics(sys, 0, an.tend, an.stepSize, an.outputSteps);
end
fprintf('time for analysis: %f\n', toc);

%% Postprocessing of simulation data...
vis = input('Animate? Y/N [Y]:', 's');
if isempty(vis)
    vis = 'Y';
end

if strcmp(vis, 'Y')
    dlmwrite('foo.out', [data.t ; data.q]');
    visualize(model_file, 'foo.out');
end
% TODO