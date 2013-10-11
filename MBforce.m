classdef MBforce < matlab.mixin.Heterogeneous
    %MBforce Base class for force elements
    %   Detailed explanation goes here
    
    properties
        id
        bodyI
        bodyJ = [];
    end
    
    methods
        function obj = MBforce(id, bodyI, bodyJ)
            if nargin > 0
                obj.id = id;
                obj.bodyI = bodyI;
                obj.bodyJ = bodyJ;
            end
        end
        
        function print(obj)
            fprintf('ForceElement %i    type: %s\n', obj.id, obj.getType);
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

