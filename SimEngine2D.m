function [sys, data] = SimEngine2D(model_file, analysis_file)
% SimEngine2D - dynamics and kinematics analysis of planar mechanisms.

%% Create the MBsystem using the definition in the ADM file.
sys = MBsys(model_file);

%% Parse the ACF file to extract analysis information.
an = loadjson(analysis_file);

if strcmpi(an.simulation, 'kinematics')
    data = kinematics(sys, 0, an.tend, an.stepSize, an.outputSteps);
elseif strcmpi(an.simulation, 'dynamics')
    data = dynamics(sys, 0, an.tend, an.stepSize, an.outputSteps);
end

%% Postprocessing of simulation data...
% TODO