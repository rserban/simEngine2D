function data = analyze(filenameACF, sys)

fid = fopen(filenameACF);
C = textscan(fid, '%s%s', 'CommentStyle', '%');
fclose(fid);

keys = C{1};
values = C{2};

for it = 1:length(keys)
    if strcmp(keys{it}, 'simulation:')
        sim = values{it};
    elseif strcmp(keys(it), 'tend:')
        tend = str2double(values{it});
    elseif strcmp(keys{it}, 'stepSize:')
        dt = str2double(values{it});
    elseif strcmp(keys{it}, 'outputSteps:')
        nOutSteps = str2double(values{it});
    else
        warning('Unknown token');
    end
end

if strcmp(sim, 'kinematics')
    data = kinematics(sys, 0, tend, dt, nOutSteps);
end