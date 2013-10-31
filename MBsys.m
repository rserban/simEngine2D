classdef MBsys < handle
    %MBSYS This class encapsulates a multibody system.
    %   Detailed explanation goes here
    
    properties
        name = 'model';   % name of this model
        nB = 0;           % number of bodies
        nC = 0;           % number of constraint elements (joints)
        nF = 0;           % number of force elements
        n = 0;            % number of generalized states
        m = 0;            % number of constraint equations
        bodies            % list of bodies
        constraints       % list of constraint elements (joints)
        forces            % list of force elements
        
        g = [0; -9.81];   % gravitatinal acceleration
    end
    
    methods
        function obj = MBsys(filename)
            % Construct an MBsys from the specification in a JSON file.
            model = loadjson(filename);
            
            if isfield(model, 'name')
                obj.name = model.name;
            else
                obj.name = 'model';
            end
            if isfield(model, 'gravity')
                obj.g = model.gravity';
            end
            
            % Create the bodies...
            if isfield(model, 'bodies')
                if iscell(model.bodies)
                    blist = model.bodies;
                else
                    blist = num2cell(model.bodies);
                end
                obj.nB = length(blist);
                obj.bodies = MBbody.empty(obj.nB,0);
            end
            
            obj.n = 0;
            for iB = 1:obj.nB
                obj.bodies(iB) = MBbody(blist{iB});
                obj.bodies(iB).setStartIndex(obj.n+1);
                obj.n = obj.n + 3;
            end
            
            % Create the constraints...
            if isfield(model, 'constraints')
                if iscell(model.constraints)
                    clist = model.constraints;
                else
                    clist = num2cell(model.constraints);
                end
                obj.nC = length(clist);
                obj.constraints = MBconstraint.empty(obj.nC,0);
            end
            
            obj.m = 0;
            for iC = 1:obj.nC
                if strcmp(clist{iC}.type, 'AbsoluteX')
                    obj.constraints(iC) = MBabsX(clist{iC});
                elseif strcmp(clist{iC}.type, 'AbsoluteY')
                    obj.constraints(iC) = MBabsY(clist{iC});
                elseif strcmp(clist{iC}.type, 'AbsoluteAngle')
                    obj.constraints(iC) = MBabsAngle(clist{iC});
                elseif strcmp(clist{iC}.type, 'RevoluteJoint')
                    obj.constraints(iC) = MBrevJoint(clist{iC});
                end
                obj.constraints(iC).start = obj.m+1;
                obj.m = obj.m + obj.constraints(iC).neq;
            end
            
            % Create the force elements...
            if isfield(model, 'forces')
                if iscell(model.forces)
                    flist = model.forces;
                else
                    flist = num2cell(model.forces);
                end
                obj.nF = length(flist);
                obj.forces = MBforce.empty(obj.nF,0);
            end
            
            for iF = 1:obj.nF
                if strcmp(flist{iF}.type, 'PointForce')
                    obj.forces(iF) = MBpointForce(flist{iF});
                elseif strcmp(flist{iF}.type, 'Torque')
                    obj.forces(iF) = MBtorque(flist{iF});
                elseif strcmp(flist{iF}.type, 'RSDA1')
                    obj.forces(iF) = MBrsda1(flist{iF});
                elseif strcmp(flist{iF}.type, 'RSDA2')
                    obj.forces(iF) = MBrsda2(flist{iF});
                end
            end
        end
        
        function [q0, qd0] = setIC(obj)
            q0 = zeros(obj.n, 1);
            qd0 = zeros(obj.n, 1);
            
            for iB = 1:obj.nB
                body = obj.bodies(iB);
                q0(body.start : body.start+2) = body.q0;
                qd0(body.start : body.start+2) = body.qd0;
            end
        end
        
        function M = getM(obj)
            M = zeros(obj.n, obj.n);
            
            for iB = 1:obj.nB
                body = obj.bodies(iB);
                M(body.start,body.start) = body.mass;
                M(body.start+1,body.start+1) = body.mass;
                M(body.start+2,body.start+2) = body.jbar;
            end
        end
        
        function Q = evalQ(obj, t, q, qd)
            Q = zeros(obj.n,1);
            
            % Apply gravity to all bodies in the system.
            for iB = 1:obj.nB
                body = obj.bodies(iB);
                Q(body.start:body.start+1) = body.mass * obj.g;
            end
            
            % Allow all force elements to add their contributions.
            for iF = 1:obj.nF
                force = obj.forces(iF);
                bodyI = obj.bodies(force.bodyI);
                i_range = (bodyI.start : bodyI.start+2);
                if isempty(force.bodyJ)
                    Qi = force.eval(t, q(i_range), qd(i_range));
                    Q(i_range) = Q(i_range) + Qi;
                else
                    bodyJ = obj.bodies(force.bodyJ);
                    j_range = (bodyJ.start : bodyJ.start + 2);
                    [Qi, Qj] = force.eval(t, q(i_range), qd(i_range), q(j_range), qd(j_range));
                    Q(i_range) = Q(i_range) + Qi;
                    Q(j_range) = Q(j_range) + Qj;
                end
            end
        end
        
        function Phi = evalPhi(obj, t, q)
            Phi = zeros(obj.m, 1);
            
            for iC = 1:obj.nC
                cnstr = obj.constraints(iC);
                eqs_range = (cnstr.start : cnstr.start+cnstr.neq-1);
                bodyI = obj.bodies(cnstr.bodyI);
                qI_range = (bodyI.start : bodyI.start+2);
                if isempty(cnstr.bodyJ)
                    [Phi(eqs_range),~,~,~] = cnstr.eval(t, q(qI_range), [], [1,0,0,0]);
                else
                    bodyJ = obj.bodies(cnstr.bodyJ);
                    qJ_range = (bodyJ.start : bodyJ.start+2);
                    [Phi(eqs_range),~,~,~] = cnstr.eval(t, q(qI_range), q(qJ_range), [], [], [1,0,0,0]);
                end
            end
        end
        
        function Phi_q = evalPhi_q(obj, t, q)
            Phi_q = zeros(obj.m, obj.n);
            
            for iC = 1:obj.nC
                cnstr = obj.constraints(iC);
                eqs_range = (cnstr.start : cnstr.start+cnstr.neq-1);
                bodyI = obj.bodies(cnstr.bodyI);
                qI_range = (bodyI.start : bodyI.start+2);
                if isempty(cnstr.bodyJ)
                    [~,Phi_q(eqs_range,qI_range),~,~] = cnstr.eval(t, q(qI_range), [], [0,1,0,0]);
                else
                    bodyJ = obj.bodies(cnstr.bodyJ);
                    qJ_range = (bodyJ.start : bodyJ.start+2);
                    [~,cPhi_q,~,~] = cnstr.eval(t, q(qI_range), q(qJ_range), [], [], [0,1,0,0]);
                    Phi_q(eqs_range,qI_range) = cPhi_q(:,1:3);
                    Phi_q(eqs_range,qJ_range) = cPhi_q(:,4:6);
                end
            end
        end
        
        function Nu = evalNu(obj, t, q)
            Nu = zeros(obj.m, 1);
            
            for iC = 1:obj.nC
                cnstr = obj.constraints(iC);
                eqs_range = (cnstr.start : cnstr.start+cnstr.neq-1);
                bodyI = obj.bodies(cnstr.bodyI);
                qI_range = (bodyI.start : bodyI.start+2);
                if isempty(cnstr.bodyJ)
                    [~,~,Nu(eqs_range),~] = cnstr.eval(t, q(qI_range), [], [0,0,1,0]);
                else
                    bodyJ = obj.bodies(cnstr.bodyJ);
                    qJ_range = (bodyJ.start : bodyJ.start+2);
                    [~,~,Nu(eqs_range),~] = cnstr.eval(t, q(qI_range), q(qJ_range), [], [], [0,0,1,0]);
                end
            end
        end
        
        function Gamma = evalGamma(obj, t, q, qd)
            Gamma = zeros(obj.m, 1);
            
            for iC = 1:obj.nC
                cnstr = obj.constraints(iC);
                eqs_range = (cnstr.start : cnstr.start+cnstr.neq-1);
                bodyI = obj.bodies(cnstr.bodyI);
                qI_range = (bodyI.start : bodyI.start+2);
                if isempty(cnstr.bodyJ)
                    [~,~,~,Gamma(eqs_range)] = cnstr.eval(t, q(qI_range), qd(qI_range), [0,0,0,1]);
                else
                    bodyJ = obj.bodies(cnstr.bodyJ);
                    qJ_range = (bodyJ.start : bodyJ.start+2);
                    [~,~,~,Gamma(eqs_range)] = cnstr.eval(t, q(qI_range), q(qJ_range), qd(qI_range), qd(qJ_range), [0,0,0,1]);
                end
            end
        end
        
        function print(obj)
            fprintf('\nModel name: %s\n', obj.name);
            fprintf('Gravity: (%g  %g)\n', obj.g);
            
            fprintf('\nNumber of bodies: %i\n', obj.nB);
            for iB = 1:obj.nB
                obj.bodies(iB).print;
            end
            
            fprintf('\nNumber of constraints: %i\n', obj.nC);
            for iC = 1:obj.nC
                obj.constraints(iC).print;
            end
            
            fprintf('\nNumber of force elements: %i\n', obj.nF);
            for iF = 1:obj.nF
                obj.forces(iF).print;
            end
            
            fprintf('\nNumber of generalized states:   %i\n', obj.n);
            fprintf('Number of constraint equations: %i\n', obj.m);
        end
        
        function plotBody(obj, data, bodyId)
            % Plot time evolution for the generalized coordinates and their
            % derivatives for the specified body, using the provided time
            % history for generalized coordinates and derivatives.
            figName = sprintf('[%s]  Body #%i', obj.name, bodyId);
            hf = figure;
            set(hf, 'position', [500, 100, 560, 640]);
            set(hf, 'name', figName, 'numbertitle', 'off');
            
            body = obj.bodies(bodyId);
            range = (body.start:body.start+2);
            
            subplot(3,1,1), hold on, box on, grid on, title('Generalized coordinates')
            plot(data.t, data.q(range(1),:), 'r');
            plot(data.t, data.q(range(2),:), 'g');
            plot(data.t, data.q(range(3),:), 'b');
            legend('x', 'y', '\phi');
            
            subplot(3,1,2), hold on, box on, grid on, title('Generalized velocities')
            plot(data.t, data.qd(range(1),:), 'r');
            plot(data.t, data.qd(range(2),:), 'g');
            plot(data.t, data.qd(range(3),:), 'b');
            legend('x', 'y', '\phi');
            
            subplot(3,1,3), hold on, box on, grid on, title('Generalized accelerations')
            plot(data.t, data.qdd(range(1),:), 'r');
            plot(data.t, data.qdd(range(2),:), 'g');
            plot(data.t, data.qdd(range(3),:), 'b');
            legend('x', 'y', '\phi');
            
        end
        
        function plotPoint(obj, data, bodyId, P)
            % Plot time evolution for the specified point on the specified
            % body, using the provided history for generalized coordinates
            % and derivatives.
            
            body = obj.bodies(bodyId);
            range = (body.start:body.start+2);
            
            q = data.q(range,:);
            qd = data.qd(range,:);
            qdd = data.qdd(range,:);
            c = cos(q(3,:));
            s = sin(q(3,:));
            r = [q(1,:) + P(1)*c - P(2)*s;
                q(2,:) - P(1)*s + P(2)*c];
            rd = [qd(1,:) - P(1)*s.*qd(3,:) - P(2)*c.*qd(3,:);
                qd(2,:) - P(1)*c.*qd(3,:) - P(2)*s.*qd(3,:)];
            rdd = [qdd(1,:) - P(1)*s.*qdd(3,:) - P(1)*c.*qd(3,:).^2 - P(2)*c.*qdd(3,:) + P(2)*s.*qd(3,:).^2;
                qdd(2,:) - P(1)*c.*qdd(3,:) + P(1)*s.*qd(3,:).^2 - P(2)*s.*qdd(3,:) - P(2)*c.*qd(3,:).^2];
            
            % Plot time evolution of position, velocity, acceleration.
            figName = sprintf('[%s]  Point (%g,%g) on Body #%i', obj.name, P(1), P(2), bodyId);
            hf1 = figure;
            set(hf1, 'position', [500, 100, 560, 640]);
            set(hf1, 'name', figName, 'numbertitle', 'off');
            
            subplot(3,1,1), hold on, box on, grid on, title('Position')
            plot(data.t, r(1,:), 'r');
            plot(data.t, r(2,:), 'g');
            legend('x', 'y');
            
            subplot(3,1,2), hold on, box on, grid on, title('Velocity')
            plot(data.t, rd(1,:), 'r');
            plot(data.t, rd(2,:), 'g');
            legend('x', 'y');
            
            subplot(3,1,3), hold on, box on, grid on, title('Acceleration')
            plot(data.t, rdd(1,:), 'r');
            plot(data.t, rdd(2,:), 'g');
            legend('x', 'y');
            
            
            % Plot trajectory and phase-plot.
            hf2 = figure;
            set(hf2, 'position', [550, 100, 560, 640]);
            set(hf2, 'name', figName, 'numbertitle', 'off');
            
            subplot(2,1,1), hold on, box on, grid on, title('Trajectory')
            plot(r(1,:), r(2,:), 'k')
            plot(r(1,1), r(2,1), 'ko')
            xlabel('x')
            ylabel('y')
            axis square
            
            subplot(2,1,2), hold on, box on, grid on, title('Phase-plot')
            plot(r(1,:), rd(1,:), 'r')
            plot(r(2,:), rd(2,:), 'g')
            legend('x', 'y');
            plot(r(1,1), rd(1,1), 'ro')
            plot(r(2,1), rd(2,1), 'go')
            xlabel('position');
            ylabel('velocity');
            
        end
        
        function [F,T] = plotReaction(obj, data, cnstrId, bodyId, P)
            % Sanity check
            cnstr = obj.constraints(cnstrId);
            
            if cnstr.bodyI == bodyId
                which = 1;
            elseif cnstr.bodyJ == bodyId
                which = 2;
            else
                warning('Inconsistent constraint/body pair.');
                return;
            end
            
            if ~isfield(data, 'lam')
                warning('The provided data structure does not contain Lagrange multipliers.');
                return;
            end
            
            % Extract states for bodyI and bodyJ (as needed).
            bodyI = obj.bodies(cnstr.bodyI);
            i_range = (bodyI.start:bodyI.start+2);
            qi = data.q(i_range,:);
            
            if ~isempty(cnstr.bodyJ)
                bodyJ = obj.bodies(cnstr.bodyJ);
                j_range = (bodyJ.start:bodyJ.start+2);
                qj = data.q(j_range,:);
            end
            
            % Extract Lagrange multipliers for this constraint.
            l_range = (cnstr.start:cnstr.start+cnstr.neq-1);
            lam = data.lam(l_range,:);
            
            % Preallocate space for the reaction force and torque.
            nt = length(data.t);
            F = zeros(2, nt);
            T = zeros(1, nt);
            
            % Loop over all output times in the specified data and
            % calculate reaction forces.
            if isempty(cnstr.bodyJ)
                for it = 1:nt
                    [F(:,it), T(:,it)] = cnstr.calcReaction(data.t(it), 1, P, qi(:,it), [], lam(:,it));
                end
            else
                for it = 1:nt
                    [F(:,it), T(:,it)] = cnstr.calcReaction(data.t(it), which, P, qi(:,it), qj(:,it), lam(:,it));
                end
            end
            
            % Plot force and torque
            figName = sprintf('[%s]  Reactions at (%g,%g) on Body #%i due to Constraint #%i', obj.name, P(1), P(2), bodyId, cnstrId);
            hf1 = figure;
            set(hf1, 'position', [500, 100, 560, 640]);
            set(hf1, 'name', figName, 'numbertitle', 'off');
            
            subplot(2,1,1), hold on, box on, grid on, title('Reaction forces')
            plot(data.t, F(1,:), 'r');
            plot(data.t, F(2,:), 'g');
            legend('F^R_x', 'F^R_y');
            
            subplot(2,1,2), hold on, box on, grid on, title('Reaction torque')
            plot(data.t, T, 'b');
            legend('T^R')
        end
        
        function [KE, PE] = plotEnergy(obj, data)
            % Preallocate space
            nt = length(data.t);
            KE = zeros(1,nt);
            PE = zeros(1,nt);
            
            % Calculate kinetic energy.
            M = obj.getM;
            for it = 1:nt
                KE(it) = 0.5 * data.qd(:,it)' * M * data.qd(:,it);
            end
            
            % Calculate potential energy.
            J_locs = (1:obj.nB) * 3;
            M(J_locs,J_locs) = 0;
            g_ext = repmat([obj.g;0]', 1, obj.nB);
            for it = 1:nt
                PE(it) = -g_ext * M * data.q(:,it);
            end
            
            % Plots.
            figName = sprintf('[%s]  System energy', obj.name);
            hf1 = figure;
            set(hf1, 'position', [500, 100, 400, 400]);
            set(hf1, 'name', figName, 'numbertitle', 'off');
            
            hold on, box on, grid on
            plot(data.t, KE, 'b');
            plot(data.t, PE, 'r');
            plot(data.t, KE+PE, 'g');
            legend('Kinetic', 'Potential', 'Total');
        end
    end
    
end

