classdef Integrator < handle & TestEnv
    %   INTEGRATOR 
    
    properties
        dyn;   % handle for the classDyn object providing the right hand side of the ODE
        solver; %function handle to the solver
    end
    methods
        %constructor
        function fE = Integrator(varargin)
            if(nargin == 0)
                global TEST;
                
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
                % constructor based on two input values
                % a classDyn element and a classOCPparam element
            elseif (nargin == 1)
                if (isa(varargin{1},'Dyn'))
                    fE.dyn = varargin{1};
                else
                    error('wrong class type for robot dynamics');
                end
             else
                error('wrong number of inputs');
            end
        end
        
        function ode(obj)
            
        end
        
        function odeD(obj)
            
        end
        
        function odeDD(obj)
            
        end
    end
end