classdef Constraint < handle
% classConstraints providing discretized ODE constraint using forward euler

    properties
        vec;    % optimization vector combining all control and state values
        dode;   % handle for the classForwEuler element providing the discretization of the ode
        
        dyn;    % Dynamik
    end
    
    methods
        %constructor
        function cC = Constraint(varargin)
            % constructor based on two input values
            % a classForwEuler element and a classOCPparam element
            if (nargin == 2)
                if (isa(varargin{2},'ForwEuler'))
                    cC.dode = varargin{2};
                    cC.dyn = cC.dode.dyn; % Initializierung von Dynamik
                else
                    error('wrong class type for discretized ode');
                end
            else
                error('wrong number of inputs');     
            end
        end
        
        %set methods
        function set.vec(obj,vec)
            % set new input vector
            obj.vec = vec;
        end
        
        % other functions
        function [ineq_con,eq_con,ineq_conD,eq_conD] = constr(obj)
            % provide the equality and inequality constraints (along
            % with their jacobians by calling get_eq_con and 
            % get_eq_conD
            ineq_con    = [];
            eq_con      = obj.get_eq_con();
            ineq_conD   = [];
            eq_conD     = obj.get_eq_conD();
        end
        
        function eq_con = get_eq_con(obj)
            % the equality constraint of the ocp
            % combine the discretized ode with the boundary conditions
            xbc         = obj.dyn.enviroment.xbc;
            state_mat   = obj.dyn.state;

            eq_con      = [obj.dode.h(); state_mat(:,1) - xbc(:,1); ...
                                state_mat(:,end) - xbc(:,end)];
        end
        
        function eq_conD = get_eq_conD(obj)
            % the Jacobian of the equality contraints of the ocp
            [n_int, n_state, n_contr] = getParams(obj);

            srow        = 1:2*n_state;
            scol        = [1:n_state,...
                n_int*(n_state+n_contr)+1:n_int*(n_state+n_contr)+n_state];
            sval        = ones(1,2*n_state);

            eq_conD     = [obj.dode.hD(); sparse(srow,scol,sval,...
                                2*n_state,(n_int+1)*(n_state+n_contr))]';

        end
        
        function eq_conDD = get_eq_conDD(varargin)
            % the Hessian of the equality constraints of the ocp.
            % Returns a cell array of Hessians if only the object 
            % is provided. If the object and the Langrange multipliers are
            % provided, the function returns the Hessian of the Lagrangian
            if (nargin == 1)
                obj         = varargin{1};
                [n_int, n_state, n_contr] = getParams(obj);
            
                eq_conDD = [obj.dode.hDD(); cell(2*n_state,1)];
                for i=0:2*n_state-1
                    eq_conDD{end-i} = sparse((n_int+1)*(n_state+n_contr),...
                        (n_int+1)*(n_state+n_contr));
                end
            elseif (nargin == 2)
                lambda      = varargin{2}; 
                H           = varargin{1}.dode.hDD();
                
                if (~isempty(lambda))
                    eq_conDD    = {lambda(1)*H{1}};
                    for i=2:length(lambda)
                        eq_conDD{1} = eq_conDD{1} + lambda(i)*H{i}; 
                    end
                else
                    eq_conDD = cell(1,1);
                end
            else 
                eq_conDD = cell(1,1);
            end
        end
        
        function [n_int, n_state, n_contr] = getParams(obj)
            n_int       = obj.dyn.environment.n_intervals;
            n_state     = obj.dyn.robot.n_state;
            n_contr     = obj.dyn.robot.n_contr;
        end
        
        % general type get functions -> interace for testing
        function f = get_func(obj)
            % interfacing get_eq_con
            f = get_eq_con(obj);
        end
        
        function g = get_jac(obj)
            % interfacing get_eq_conD
            g = get_eq_conD(obj);
        end
        
        function H = get_hess(varargin)
            % interfacing get_eq_conDD
            if (nargin == 1)
                H = get_eq_conDD(varargin{1});
            elseif (nargin == 2)
                H = get_eq_conDD(varargin{1},varargin{2});
            else
                H = [];
            end
        end
 
    end
end