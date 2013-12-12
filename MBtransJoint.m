classdef MBtransJoint < MBconstraint
    
    properties
        sPI
        sPJ
        vPI
        vPJ
    end
    
    methods
        function obj = MBtransJoint(varargin)
            data = varargin{1};
            obj = obj@MBconstraint(data.id, 2, data.body1, data.body2, data.fun);
            obj.sPI = reshape(data.sP1, 2, 1);
            obj.sPJ = reshape(data.sP2, 2, 1);
            obj.vPI = reshape(data.vP1, 2, 1);
            obj.vPJ = reshape(data.vP2, 2, 1);
            obj.vPI = obj.vPI / norm(obj.vPI);
            obj.vPJ = obj.vPJ / norm(obj.vPJ);
        end

  
        function print(obj)
            print@MBconstraint(obj);
            fprintf('   sPI:  (%g  %g)\n', obj.sPI(1), obj.sPI(2));
            fprintf('   sPJ:  (%g  %g)\n', obj.sPJ(1), obj.sPJ(2));
            fprintf('  ... TODO ...\n');
        end
        
        
        function [Phi, Phi_q, Nu, Gamma] = eval(obj, t, qi, qj, qdi, qdj, flags)
            ri = qi(1:2);  phii = qi(3);
            rj = qj(1:2);  phij = qj(3);
            
            [Ai, Bi] = MBbody.rotMat(phii);
            [Aj, Bj] = MBbody.rotMat(phij);
            R = [0, -1; 1, 0];

            riP = ri + Ai*obj.sPI;
            rjP = rj + Aj*obj.sPJ;
            dij = rjP - riP;

            vi = Ai * obj.vPI;
            vj = Aj * obj.vPJ;

            vi_perp = R * vi;
            
            Phi = [];
            Phi_q = [];
            Nu = [];
            Gamma = [];
            
            if flags(1)
                Phi = [ vi_perp' * dij ; vi_perp' * vj ];
            end
            if flags(2)
                Phi_qi = [ -obj.vPI' * Bi'  , -obj.vPI' * Ai' * (rj - ri + Aj * obj.sPJ)
                              zeros(1,2)    , -obj.vPI' * Ai' * Aj * obj.vPJ             ];
                Phi_qj = [  obj.vPI' * Bi'  , obj.vPI' * Ai' * Aj * obj.sPJ 
                              zeros(1,2)    , obj.vPI' * Ai' * Aj * obj.vPJ  ];
                Phi_q = [Phi_qi , Phi_qj];
            end
            if flags(3)
                Nu = [0; 0];
            end
            if flags(4)
                phidi = qdi(3);
                phidj = qdj(3);
                %
                %   TO DO !!!!!!!
                %
                Gamma = [0;0];
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
            type = 'TranslationalJoint';
        end
    end
    
end