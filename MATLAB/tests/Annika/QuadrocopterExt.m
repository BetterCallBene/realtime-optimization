classdef QuadrocopterExt < BasisQDyn
    %QUADROCOPTEREXT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function res = dot(obj,ind) % Rechte Seite der ODE aus (Quadrocoptermodell.pdf: (1.46))
            % compute the right hand side of the ode for a given time
            % instance ind
            res = BasisQDyn@dot(ind);
        end
    end
    
end

