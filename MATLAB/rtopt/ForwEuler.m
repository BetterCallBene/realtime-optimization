classdef ForwEuler < Solver
    %ForwEuler Integriert den Einschritt ForwEuler in das Realtimeprojekt 
    
    properties
    end
    
    methods
        function FE = ForwEuler(varargin)
            FE@Solver(nargin, varargin); %Bug von Matlab
        end
        
        function y = integrate(obj, func, meshGrid, y0,yp0)
            y = y0 + obj.h .* feval(func, meshGrid, y0);
        end
        
    end
    
    
    
    
end

