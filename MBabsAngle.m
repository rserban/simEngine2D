classdef MBabsAngle < MBconstraint
    
    properties
    end
    
    methods
        function obj = MBabsAngle(varargin)
            data = varargin{1};
            obj = obj@MBconstraint(data.id, 1, data.body1, [], data.fun);
        end
        
        function print(obj)
            print@MBconstraint(obj);
        end
        
        function [Phi, Phi_q, Nu, Gamma] = eval(obj, t, qi, qdi, flags)
            phi = qi(3);
         
            Phi = [];
            Phi_q = [];
            Nu = [];
            Gamma = [];

            if flags(1)
                Phi = phi - obj.fun(t);
            end
            if flags(2)
                Phi_q = [0,  0,  1];
            end
            if flags(3)
                Nu = obj.funD(t);
            end
            if flags(4)
                Gamma = obj.funDD(t);
            end
        end
        
    end
    
    methods(Static)
        function type = getType()
            type = 'AbsoluteAngle';
        end
    end
    
end