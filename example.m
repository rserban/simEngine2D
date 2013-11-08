%% Create the MBsystem using the definition in the ADM file.
sys = MBsys('Models/rodDYN.adm');

%% Parse the ACF file to extract analysis information.
an = loadjson('Models/dynamics.acf');

if strcmpi(an.simulation, 'kinematics')
    data = kinematics(sys, 0, an.tend, an.stepSize, an.outputSteps);
elseif strcmpi(an.simulation, 'dynamics')
    data = dynamics(sys, 0, an.tend, an.stepSize, an.outputSteps);
end

%% Plot results
sys.plotBody(data,1);
pause
sys.plotReaction(data,2,1,[0;-2]);
pause
sys.plotReaction(data,3,1,[0;0]);
pause
sys.plotEnergy(data);
pause

%% Postprocessing of simulation data...
vis = input('Animate? Y/N [Y]:', 's');
if isempty(vis)
    vis = 'Y';
end

if strcmp(vis, 'Y')
    dlmwrite('foo.out', [data.t ; data.q]');
    visualize('Models/rodDYN.adm', 'foo.out');
end
% TODO