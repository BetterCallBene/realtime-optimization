classdef ode45M < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        function M45 = ode45M(varargin)
            M45@Solver(nargin, varargin); %Bug in Matlab
        end
        
        function y = integrate(obj, func, meshGrid, y0, yp0)
            opts_ = obj.opts;
            [t, y] = ode45(func, meshGrid, y0, opts_);
            y = y(end, :)';
        end
        
        
    end
end

