classdef MBrevJoint < MBconstraint
    
    properties
        sI
        sJ
    end
    
    methods
        function obj = MBrevJoint(varargin)
            data = varargin{1};
            obj = obj@MBconstraint(data.id, 2, data.body1, data.body2, data.fun);
            obj.sI = reshape(data.sP1, 2, 1);
            obj.sJ = reshape(data.sP2, 2, 1);
        end

  
        function print(obj)
            print@MBconstraint(obj);
            fprintf('   sI:  (%g  %g)\n', obj.sI(1), obj.sI(2));
            fprintf('   sJ:  (%g  %g)\n', obj.sJ(1), obj.sJ(2));
        end
        
        
        function [Phi, Phi_q, Nu, Gamma] = eval(obj, t, qi, qj, qdi, qdj, flags)
            ri = qi(1:2);  phii = qi(3);
            rj = qj(1:2);  phij = qj(3);
            
            [Ai, Bi] = MBbody.rotMat(phii);
            [Aj, Bj] = MBbody.rotMat(phij);
            
            Phi = [];
            Phi_q = [];
            Nu = [];
            Gamma = [];
            
            if flags(1)
                Phi = ri + Ai * obj.sI - rj - Aj * obj.sJ;
            end
            if flags(2)
                Phi_q = [eye(2), Bi * obj.sI, -eye(2), -Bj * obj.sJ];
            end
            if flags(3)
                Nu = [0; 0];
            end
            if flags(4)
                phidi = qdi(3);
                phidj = qdj(3);
                Gamma = Ai * obj.sI * phidi^2 - Aj * obj.sJ * phidj^2;
            end
        end
    
        
        function [] = evalAndPrint(obj, t, qi, qj, qdi, qdj)
            [Phi, Phi_q, Nu, Gamma] = obj.eval(t, qi, qj, qdi, qdj, [1,1,1,1]);
            fprintf('Phi    = [ %10.4f ]\n         [ %10.4f ]\n', Phi(1), Phi(2));
            fprintf('Phi_qi = [ %10.4f  %10.4f  %10.4f ]\n         [ %10.4f  %10.4f  %10.4f ]\n',...
                Phi_q(1:2,1:3)');
            fprintf('Phi_qj = [ %10.4f  %10.4f  %10.4f ]\n         [ %10.4f  %10.4f  %10.4f ]\n',...
                Phi_q(1:2,4:6)');
            fprintf('Nu     = [ %10.4f ]\n         [ %10.4f ]\n', Nu(1), Nu(2));
            fprintf('Gamma  = [ %10.4f ]\n         [ %10.4f ]\n', Gamma(1), Gamma(2));            
        end

    end
        
    methods(Static)
        function type = getType()
            type = 'RevoluteJoint';
        end
    end
    
end