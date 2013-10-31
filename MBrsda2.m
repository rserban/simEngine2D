classdef MBrsda2 < MBforce
    %MBrsda2 RSDA force element between two bodies
    %   It is assumed that an appropriate RevoluteJoint constraint
    %   is provided to pin the two bodies.
    
    properties
        theta0       % RSDA free angle
        k            % RSDA spring constant
        c            % RSDA damping coefficient
        funH;        % RSDA actuator function h(t)
        matlab_file  % name of a matlab file defining h(t,theta,thetaD)
    end
    
    methods
        function obj = MBrsda2(varargin)
            data = varargin{1};
            
            obj = obj@MBforce(data.id, data.body1, data.body2);
            
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
                
        function [Qi, Qj] = eval(obj, t, qi, qdi, qj, qdj)
            n = obj.k * (qj(3) - qi(3) - obj.theta0) + obj.c * (qdj(3) - qdi(3)) + obj.funH(t);
            
            Qi = [0; 0;  n];
            Qj = [0; 0; -n];
        end
        
    end
    
    methods(Static)
        function type = getType()
            type = 'RSDA2';
        end
    end
    
end