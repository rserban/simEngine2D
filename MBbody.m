classdef MBbody < handle
    %MBbody This class encapulates a body in 2D.
    %   Detailed explanation goes here
    
    properties
        id                % id of this body in the system
        start = 0;        % start index for this body's states in global arrays
        mass              % mass
        jbar              % Jzz (expressed in local frame)
        q0                % initial generalized states (r and phi)
        qd0               % initial generalized velocities (rd and phid)
    end
    
    methods
        function body = MBbody(varargin)
            if nargin == 1
                data = varargin{1};
                body.id = data.id;
                body.mass = data.mass;
                body.jbar = data.jbar;
                body.q0 = reshape(data.q0, 3, 1);
                body.qd0 = reshape(data.qd0, 3, 1);
            end
        end
        
        function obj = setStartIndex(obj, start)
            obj.start = start;
        end
        
        function print(obj)
            fprintf('Body %i\n', obj.id);
            fprintf('   start: %i\n', obj.start);
            fprintf('   mass: %g   jbar: %g\n', obj.mass, obj.jbar);
            fprintf('   initial gen. states:   %g  %g  %g\n', obj.q0(1), obj.q0(2), obj.q0(3));
            fprintf('   init gen. velocities:  %g  %g  %g\n', obj.qd0(1), obj.qd0(2), obj.qd0(3));
        end
    end
    
    methods(Static)
        function [A, varargout] = rotMat(angle)
            c = cos(angle);
            s = sin(angle);
            
            A = [c -s; s c];
            if (nargout == 2)
                varargout{1} = [-s -c; c -s];
            end
        end
    end
end

