classdef Sensitivity < handle & TestEnv
    %   INTEGRATOR 
    
    properties
        dyn;   % handle for the classDyn object providing the right hand side of the ODE
        solver; %function handle to the solver
    end
    methods(Access=private)
        
    end
    methods
        %constructor
        function Se = Sensitivity(varargin)
            if(nargin == 0)
                global TEST;
                
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
                % constructor based on two input values
                % a classDyn element and a classOCPparam element
            elseif (nargin == 2)
                if (isa(varargin{1},'Dyn'))
                    Se.dyn = varargin{1};
                else
                    error('wrong class type for robot dynamics');
                end
                if (isa(varargin{2},'Solver'))
                    Se.solver = varargin{2};
                else
                    error('wrong class type for robot dynamics');
                end
             else
                error('wrong number of inputs');
            end
        end
        
        function res = ode(obj, h, y)
            solver(  
        end
        
    end
end