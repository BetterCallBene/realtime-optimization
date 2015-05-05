classdef classEnvConstant < classEnv
    %CLASSENVCONSTANT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Access=public)
        g = 9.81;
    end
    
    methods
        function ret = get.g(obj)
            ret = obj.g;
        end
    end
    
end

