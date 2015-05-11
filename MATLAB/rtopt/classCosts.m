classdef classCosts < handle
% classCost providing cost function based on control only
%   Approximate the integral over the squared two-norm of the controls 
%   by a piecewise constant function using the discretized control values,
%   i.e., the cost value is a sum of the the squared control values
%   multiplied by the interval widths.

    properties
        vec;    % optimization vector combining all control and state values
        param;  % handle for a classOCPparam element providing the parameters describing the ocp
    end
    
    properties  (Dependent)
        contr;  % matrix of control values extracted from vec
    end
    
    methods
        %constructor
        function cC = classCosts(varargin)
            % constructor based on one input value: a classOCPparam object
            if (nargin == 1)          
                if (isa(varargin{1},'classOCPparam'))
                    cC.param = varargin{1};
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
            obj.vec = vec;
        end
        
        %get methods        
        function val = get.contr(obj)
           % get the current control values
           % (if not yet stored, extract them from vec)  
           if ((~isempty(obj.param)))
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
        function c_val = get_cost(obj)
            % compute the cost value
            control = obj.contr;
            mesh    = obj.param.mesh;
            c_val   = .5*sum((control(:,1:end-1).^2)*mesh');
        end
        
        function cD_val = get_costD(obj)
            % compute the Jacobian of the cost value
            n_tp        = obj.param.n_timepoints;
            n_var       = obj.param.n_var;
            n_contr     = obj.param.n_contr;
            control     = obj.contr;
            mesh        = obj.param.mesh;            
            
            rind        = zeros(1,(n_tp-1)*n_contr);
            vind        = zeros(1,(n_tp-1)*n_contr);
            for i = 1:n_tp-1
                rind((i-1)*n_contr+1:i*n_contr) = ...
                    (i-1)*(n_var+n_contr)+n_var + 1:...
                    i*(n_var+n_contr);
                vind((i-1)*n_contr+1:i*n_contr) = ...
                    mesh(i)*control(:,i)';
            end

            cD_val        = sparse(rind,ones(1,(n_tp-1)*n_contr),vind,...
                n_tp*(n_var+n_contr),1);          
        end
        
        function cDD_val = get_costDD(obj)
            %compute the Hessian of the cost value
            n_tp        = obj.param.n_timepoints;
            n_var       = obj.param.n_var;
            n_contr     = obj.param.n_contr;
            mesh        = obj.param.mesh;            
            
            rind        = zeros(1,(n_tp-1)*n_contr);
            vind        = zeros(1,(n_tp-1)*n_contr);
            for i = 1:n_tp-1
                rind((i-1)*n_contr+1:i*n_contr) = ...
                    (i-1)*(n_var+n_contr)+n_var + 1:...
                    i*(n_var+n_contr);
                vind((i-1)*n_contr+1:i*n_contr) = ...
                    mesh(i)*ones(1,n_contr);
            end

            cDD_val        = {sparse(rind,rind,vind,...
                n_tp*(n_var+n_contr),n_tp*(n_var+n_contr))};  
        end
        
        
        % general type get functions -> interace for testing
        function f = get_func(obj)
            f = get_cost(obj);
        end
        
        function g = get_jac(obj)
            g = get_costD(obj);
        end
        
        function H = get_hess(obj)
            H = get_costDD(obj);
        end
 
    end
end