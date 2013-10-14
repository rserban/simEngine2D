classdef MBabsX < MBconstraint
    
    properties
        sI
    end
    
    methods
        function obj = MBabsX(varargin)
            data = varargin{1};
            obj = obj@MBconstraint(data.id, 1, data.body1, [], data.fun);
            obj.sI = reshape(data.sP1, 2, 1);
        end
        
        function print(obj)
            print@MBconstraint(obj);
            fprintf('   s:  (%g  %g)\n', obj.sI(1), obj.sI(2));
        end
        
        function [Phi, Phi_q, Nu, Gamma] = eval(obj, t, qi, qdi, flags)
            x = qi(1);  phi = qi(3);
            c = cos(phi); s = sin(phi);
            
            Phi = [];
            Phi_q = [];
            Nu = [];
            Gamma = [];
            
            if flags(1)
                Phi = x + obj.sI(1)*c - obj.sI(2)*s - obj.fun(t);
            end
            if flags(2)
                Phi_q = [1,  0,  -obj.sI(1)*s - obj.sI(2)*c];
            end
            if flags(3)
                Nu = obj.funD(t);
            end
            if flags(4)
                phid = qdi(3);
                Gamma = (obj.sI(1)*c - obj.sI(2)*s) * phid^2 + obj.funDD(t);
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
            type = 'AbsoluteX';
        end
    end
    
end