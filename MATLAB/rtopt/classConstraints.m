classdef classConstraints < handle
% classConstraints providing discretized ODE constraint using forward euler

    properties
        vec;    % optimization vector combining all control and state values
        param;  % handle for a classOCPparam element providing the parameters describing the ocp
        dode;   % handle for the classForwEuler element providing the discretization of the ode
    end
    
    properties  (Dependent)
        state;  % matrix of state values extracted from vec
        contr;  % matrix of control values extracted from vec
    end
    
    methods
        %constructor
        function cC = classConstraints(varargin)
            % constructor based on two input values
            % a classForwEuler element and a classOCPparam element
            if (nargin == 2)
                if (isa(varargin{1},'classForwEuler'))
                    cC.dode = varargin{1};
                else
                    error('wrong class type for discretized ode');
                end
                
                if (isa(varargin{2},'classOCPparam'))
                    cC.param = varargin{2};
                else
                    error('wrong class type for problem parameter');
                end
            else
                error('wrong number of inputs');     
            end
        end
        
        %set methods
        function set.vec(obj,vec)
            % set new input vector 
            % and pass information on th dode 
            if ((length(obj.vec)==0)||(norm(obj.vec - vec)>1e-10))
                obj.vec = vec;
                if (~isempty(obj.dode))           
                    obj.dode.state = obj.state;
                    obj.dode.contr = obj.contr;
                end
            end
        end
        
        %get methods
        function val = get.state(obj)
           % get the current state values
           % (if not yet stored, extract them from vec)
           if ((~isempty(obj.dode))&&(~isempty(obj.param)))
                n_int       = obj.param.n_intervals;
                n_var       = obj.param.n_var;
                n_contr     = obj.param.n_contr;
                
                val        = zeros(n_var,n_int+1);

                for i = 1:n_int+1
                    val(:,i) = obj.vec((i-1)*(n_var+n_contr)+1:...
                                        (i-1)*(n_var+n_contr)+n_var);
                end
            end
        end
        
        function val = get.contr(obj)
           % get the current control values
           % (if not yet stored, extract them from vec)            
           if ((~isempty(obj.dode))&&(~isempty(obj.param)))
                n_int       = obj.param.n_intervals;
                n_var       = obj.param.n_var;
                n_contr     = obj.param.n_contr;

                val        = zeros(n_contr,n_int+1);

                for i = 1:n_int+1
                    val(:,i) = obj.vec((i-1)*(n_var+n_contr)+n_var + 1:...
                                        i*(n_var+n_contr));
                end
            end
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
            xbc         = obj.param.xbc;

            state_mat   = obj.state;

            eq_con      = [obj.dode.h(); state_mat(:,1) - xbc(:,1); ...
                                state_mat(:,end) - xbc(:,end)];
        end
        
        function eq_conD = get_eq_conD(obj)
            % the Jacobian of the equality contraints of the ocp
            n_int       = obj.param.n_intervals;
            n_var       = obj.param.n_var;
            n_contr     = obj.param.n_contr;

            srow        = 1:2*n_var;
            scol        = [1:n_var,...
                n_int*(n_var+n_contr)+1:n_int*(n_var+n_contr)+n_var];
            sval        = ones(1,2*n_var);

            eq_conD     = [obj.dode.hD(); sparse(srow,scol,sval,...
                                2*n_var,(n_int+1)*(n_var+n_contr))]';

        end
        
        function eq_conDD = get_eq_conDD(varargin)
            % the Hessian of the equality constraints of the ocp.
            % Returns a cell array of Hessians if only the object 
            % is provided. If the object and the Langrange multipliers are
            % provided, the function returns the Hessian of the Lagrangian
            if (nargin == 1)
                obj         = varargin{1};
                n_int       = obj.param.n_intervals;
                n_var       = obj.param.n_var;
                n_contr     = obj.param.n_contr;
            
                eq_conDD = [obj.dode.hDD(); cell(2*n_var,1)];
                for i=0:2*n_var-1
                    eq_conDD{end-i} = sparse((n_int+1)*(n_var+n_contr),...
                        (n_int+1)*(n_var+n_contr));
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