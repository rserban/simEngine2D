classdef MBconstraint < matlab.mixin.Heterogeneous
    %MBconstraint Base class for constraints
    %   Detailed explanation goes here
    
    properties
        id
        start = 0;
        neq
        bodyI
        bodyJ = [];
        fun
        funD
        funDD
    end
    
    methods
        function obj = MBconstraint(id, neq, bodyI, bodyJ, funStr)
            if nargin > 0
                obj.id = id;
                obj.neq = neq;
                obj.bodyI = bodyI;
                obj.bodyJ = bodyJ;
                
                % Create functions.
                % NOTE: if funStr represents a constant function (i.e. no
                % explicit dependency on 't'), eval(funStr) would result
                % in a double, and not a sym variable (since eval is an
                % overloaded function). To address this, we force the
                % result from eval(funStr) to *always* be a sym. In
                % addition, we explicitly specify the variable as 't'.
                if ~strcmp(funStr, 'NONE')
                    syms t;
                    f = sym(eval(funStr));
                    obj.fun = matlabFunction(f, 'vars', t);
                    obj.funD = matlabFunction(diff(f), 'vars', t);
                    obj.funDD = matlabFunction(diff(diff(f)), 'vars', t);
                end
            end
        end
        
        function obj = set.start(obj, start)
            obj.start = start;
        end
        
        function obj = setStartIndex(obj, start)
            obj.start = start;
        end
        
        function [F,T] = calcReaction(obj, t, which, P, qi, qj, lam)
            % Evaluate constraint Jacobian and extract Phi_r and Phi_phi.
            if isempty(obj.bodyJ)
                [~,Phi_q, ~, ~] = obj.eval(t, qi, [], [0,1,0,0]);
                [~,B] = MBbody.rotMat(qi(3));
            else
                [~,Phi_q,~,~] = obj.eval(t, qi, qj, [], [], [0,1,0,0]);
                if which == 1
                    Phi_q(:,4:6) = [];
                    [~,B] = MBbody.rotMat(qi(3));
                else
                    Phi_q(:,1:3) = [];
                    [~,B] = MBbody.rotMat(qj(3));
                end
            end
            
            F = -Phi_q(:,1:2)' * lam;
            T = -P' * B' * F - Phi_q(:,3)' * lam;
        end
        
        function print(obj)
            fprintf('Constraint %i    type: %s\n', obj.id, obj.getType);
            fprintf('   start_eq: %i  neq: %i\n', obj.start, obj.neq);
            if isempty(obj.bodyJ)
                fprintf('   body: %i\n', obj.bodyI);
            else
                fprintf('   bodyI: %i  bodyJ: %i\n', obj.bodyI, obj.bodyJ);
            end
        end
    end
    
    methods(Static)
        function type = getType()
            type = 'Generic';
        end
    end
    
end

