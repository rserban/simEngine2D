classdef MBtorque < MBforce
    
    properties
        fun;         % time function for N(t)
        matlab_file  % name of a matlab file defining the torque
    end
    
    methods
        function obj = MBtorque(varargin)
            data = varargin{1};
            obj = obj@MBforce(data.id, data.body1, []);
            syms t;
            obj.fun = matlabFunction(eval(data.fun), 'vars', t);
        end
        
        
        function print(obj)
            print@MBforce(obj);
            fprintf('   torque function:  \n');
        end
        
        function Q = eval(obj, t, qi, qdi)
            Q = [0;0;obj.fun(t)];
        end
        
    end
    
    methods(Static)
        function type = getType()
            type = 'Torque';
        end
    end
    
end