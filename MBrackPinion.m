classdef MBrackPinion < MBconstraint
    
    properties
        sPI
        sQI
        sPJ
        thetaI
        thetaJ
        vI
        vlen
        RJ
    end
    
    methods
        function obj = MBrackPinion(varargin)
            data = varargin{1};
            obj = obj@MBconstraint(data.id, 2, data.body1, data.body2, data.fun);
            obj.sPI = reshape(data.sP1, 2, 1);
            obj.sQI = reshape(data.sQ1, 2, 1);
            obj.sPJ = reshape(data.sP2, 2, 1);
            obj.thetaI = data.theta1;
            obj.thetaJ = data.theta2;
            obj.vI = obj.sQI - obj.sPI;
            obj.vlen = norm(obj.vI);
            obj.RJ = data.R2;
        end

  
        function print(obj)
            print@MBconstraint(obj);
            fprintf('   sPI:  (%g  %g)\n', obj.sPI(1), obj.sPI(2));
            fprintf('   sQI:  (%g  %g)\n', obj.sQI(1), obj.sQI(2));
            fprintf('   sPJ:  (%g  %g)\n', obj.sPJ(1), obj.sPJ(2));
            fprintf(' ... TODO ...\n');
        end
        
        
        function [Phi, Phi_q, Nu, Gamma] = eval(obj, t, qi, qj, qdi, qdj, flags)
            ri = qi(1:2);  phii = qi(3);
            rj = qj(1:2);  phij = qj(3);
            
            [Ai, Bi] = MBbody.rotMat(phii);
            [Aj, Bj] = MBbody.rotMat(phij);
            R = [0, -1; 1, 0];
            
            riP = ri + Ai*obj.sPI;
            riQ = ri + Ai*obj.sQI;
            rjP = rj + Aj*obj.sPJ;
            
            vi = Ai * obj.vI;
            vi_perp = R * vi;
            
            alphaj = phii + obj.thetaI - phij - obj.thetaJ + 3*pi/2;
            
            Phi = [];
            Phi_q = [];
            Nu = [];
            Gamma = [];
            
            if flags(1)
                Phi = [(rjP - riQ)' * vi - obj.vlen * obj.RJ * alphaj;
                       (rjP - riP)' * vi_perp - obj.vlen * obj.RJ];
            end
            if flags(2)
                Phi1_qi = [-vi' , (rjP-riQ)'*Bi*obj.vI - vi'*Bi*obj.sQI - obj.vlen*obj.RJ];
                Phi1_qj = [ vi' , vi'*Bj*obj.sPJ + obj.vlen*obj.RJ];
                Phi2_qi = [-vi_perp' , (rjP-riP)'*R*Bi*obj.vI - vi_perp'*Bi*obj.sPI];
                Phi2_qj = [ vi_perp' , vi_perp'*Bj*obj.sPJ];
                Phi_q = [Phi1_qi , Phi1_qj ; Phi2_qi , Phi2_qj];
            end
            if flags(3)
                Nu = [0; 0];
            end
            if flags(4)
                rid = qdi(1:2);   phidi = qdi(3);
                rjd = qdj(1:2);   phidj = qdj(3);
                ridP = rid + phidi * Bi * obj.sPI;
                rjdP = rjd + phidj * Bj * obj.sPJ;
                phid2_diff = phidj^2 - phidi^2;
                rd_diff = rjdP - ridP;
                Gamma1 = phid2_diff * rjP' * vi - 2 * rd_diff' * phidi * Bi * obj.vI;
                Gamma2 = phid2_diff * rjP' * vi_perp - 2 * rd_diff' * phidi * R * Bi * obj.vI;
                Gamma = [Gamma1 ; Gamma2];
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
            type = 'RackPinion';
        end
    end
    
end