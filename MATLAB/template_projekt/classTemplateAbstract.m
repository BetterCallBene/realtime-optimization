classdef classTemplateAbstract < classSTemplateAbstract
    %classTemplateAbstract Unterklasse einer abstrakten Superclass
    
    properties(Constant)
        g = 9.81;
    end
    
    methods
        function g = GetGravitation(obj)
            g = obj.g;
        end
    end
    
end

