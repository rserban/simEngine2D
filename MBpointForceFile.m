classdef MBpointForceFile < MBforce
    
    properties
        sI           % point on body where force is applied
        lrf = true;  % true if expressed in LRF, false if in GRF
        mfile;       % filename with m-function
    end
    
    methods
        function obj = MBpointForceFile(varargin)
            data = varargin{1};
            obj = obj@MBforce(data.id, data.body1, []);
            obj.sI = reshape(data.sP1, 2, 1);
            obj.mfile = data.mfile;
            if strcmp(data.frame, 'GRF')
                obj.lrf = false;
            else
                obj.lrf = true;
            end
        end
        
        
        function print(obj)
            print@MBforce(obj);
            fprintf('   sI:  (%g  %g)\n', obj.sI(1), obj.sI(2));
            if obj.lrf
                frame = 'LRF (local)';
            else
                frame = 'GRF (global)';
            end
            fprintf('   reference frame: %s\n', frame);
            fprintf('   force m-file:    %s\n', obj.mfile);
        end
        
        function Q = eval(obj, t, qi, qdi)
            % Evaluate X and Y components of the linear force.
            fstr = sprintf('%s(t, qi, qdi)', obj.mfile);
            F = eval(fstr);
            % Evaluate B matrix of the body.
            [A,B] = MBbody.rotMat(qi(3));
            % If the force is expressed in LRF, re-express it in GRF.
            if obj.lrf
                F = A * F;
            end
            % Assemble generalized force.
            Q = [F; (B*obj.sI)'*F];   
        end
        
    end
    
    methods(Static)
        function type = getType()
            type = 'PointForceFile';
        end
    end
    
end