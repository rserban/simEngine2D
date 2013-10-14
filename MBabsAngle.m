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
        
        function [] = evalAndPrint(obj, t, qi, qdi)
            [Phi, Phi_q, Nu, Gamma] = obj.eval(t, qi, qdi, [1,1,1,1]);
            fprintf('Phi    = [ %10.4f ]\n', Phi);
            fprintf('Phi_qi = [ %10.4f  %10.4f  %10.4f ]\n', Phi_q);
            fprintf('Nu     = [ %10.4f ]\n', Nu);
            fprintf('Gamma  = [ %10.4f ]\n', Gamma);  
        end
        
    end
    
    methods(Static)
        function type = getType()
            type = 'AbsoluteAngle';
        end
    end
    
end