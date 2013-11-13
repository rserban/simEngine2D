classdef MBpointForce < MBforce
    
    properties
        sI           % point on body where force is applied
        lrf = true;  % true if expressed in LRF, false if in GRF
        funX;        % time function for F_x(t)
        funY;        % time function for F_y(t)
        matlab_file  % name of a matlab file defining the force
    end
    
    methods
        function obj = MBpointForce(varargin)
            data = varargin{1};
            obj = obj@MBforce(data.id, data.body1, []);
            obj.sI = reshape(data.sP1, 2, 1);
            syms t;
            funXstr = data.funX;
            funYstr = data.funY;
            obj.funX = matlabFunction(eval(funXstr), 'vars', t);
            obj.funY = matlabFunction(eval(funYstr), 'vars', t);
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
            fprintf('   force functions: %s\n', func2str(obj.funX));
            fprintf('                    %s\n', func2str(obj.funY));
        end
        
        function Q = eval(obj, t, qi, qdi)
            % Evaluate X and Y components of the linear force.
            F = [obj.funX(t); obj.funY(t)];
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
            type = 'PointForce';
        end
    end
    
end