classdef ForwEuler < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function FE = ForwEuler(varargin)
            FE@Solver(nargin, varargin); %Bug von Matlab
        end
        
%         function res = getJDot(obj)
%             n_state = obj.dyn.robot.n_state;
%             vec = obj.vec;
%             state = vec(1:n_state);
%             u0 = vec(n_state + 1 : end);
%             res = obj.dyn.getJTilde(state, u0);
%         end
%         
%         function dy = funcToIntegrate(obj, t, varargin)
%             
%             y = varargin{1};
%             u0 = obj.u0;
%             
%             [state, M0, N0] = obj.helperCreateMatrizen(y);
%             
%             obj.vec = [state; u0];
%             FTilde = obj.dyn.FTilde(y, u0);
%             
%             kM = obj.kM(M0); 
%             kN = obj.kN(N0);
%             
%             dy = obj.helperCreateVektor(FTilde, kM, kN);
%         end
        
        function y = integrate(obj, func, meshGrid, y0,yp0)
            y = y0 + obj.h .* feval(func, meshGrid, y0);
        end
        
    end
    
    
    
    
end

