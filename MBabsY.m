classdef MBabsY < MBconstraint
    
    properties
        sI
    end
    
    methods
        function obj = MBabsY(varargin)
            data = varargin{1};
            obj = obj@MBconstraint(data.id, 1, data.body1, [], data.fun);
            obj.sI = reshape(data.sP1, 2, 1);
        end
        
        function print(obj)
            print@MBconstraint(obj);
            fprintf('   s:  (%g  %g)\n', obj.sI(1), obj.sI(2));
        end
        
        function [Phi, Phi_q, Nu, Gamma] = eval(obj, t, qi, qdi, flags)
            y = qi(2);  phi = qi(3);
            c = cos(phi); s = sin(phi);
            
            Phi = [];
            Phi_q = [];
            Nu = [];
            Gamma = [];
            
            if flags(1)
                Phi = y + obj.sI(1)*s + obj.sI(2)*c - obj.fun(t);
            end
            if flags(2)
                Phi_q = [0,  1,  obj.sI(1)*c - obj.sI(2)*s];
            end
            if flags(3)
                Nu = obj.funD(t);
            end
            if flags(4)
                phid = qdi(3);
                Gamma = (obj.sI(1)*s + obj.sI(2)*c) * phid^2 + obj.funDD(t);
            end
        end
        
    end
    
    methods(Static)
        function type = getType()
            type = 'AbsoluteY';
        end
    end
end