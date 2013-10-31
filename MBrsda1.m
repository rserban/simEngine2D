classdef MBrsda1 < MBforce
    %MBrsda1 RSDA force element between a body and ground.
    %   It is assumed that appropriate absoluteX and absoluteY
    %   constraints are provided to pin the body to ground.
    
    properties
        theta0       % RSDA free angle
        k            % RSDA spring constant
        c            % RSDA damping coefficient
        funH;        % RSDA actuator function h(t)
        matlab_file  % name of a matlab file defining h(t,theta,thetaD)
    end
    
    methods
        function obj = MBrsda1(varargin)
            data = varargin{1};
            
            obj = obj@MBforce(data.id, data.body1, []);
            
            obj.theta0 = data.theta0;
            obj.k = data.k;
            obj.c = data.c;
            
            syms t;
            funHstr = data.fun;
            obj.funH = matlabFunction(eval(funHstr), 'vars', t);
        end
        
        
        function print(obj)
            print@MBforce(obj);
            fprintf('   theta0: %g  k: %g  c: %g\n', obj.theta0, obj.k, obj.c);
        end
        
        function Q = eval(obj, t, qi, qdi)
            n = obj.k * (qi(3) - obj.theta0) + obj.c * qdi(3) + obj.funH(t);
            
            Q = [0; 0; -n];
        end
        
    end
    
    methods(Static)
        function type = getType()
            type = 'RSDA1';
        end
    end
    
end