classdef classOCPparam < handle
% classOCPparam providing parameters of OCP
    
    properties
        mesh            % time mesh (1xn) vector
        n_var           % number of state variables (int)
        n_contr         % number of control variables (int)
        dyn_param       % element of classDynParam
        xbc             % (mx2) matrix with boundary values for states
                        % at start 0 and end 1
    end
    
    properties (Dependent)
        n_intervals     % number of intervals in mesh
        n_timepoints    % number of discretization points in mesh
    end
    
    methods
        % constructor
        function cOCPp = classOCPparam(varargin)
            % constructor sets values for n_var, n_control, xbc and
            % dyn_param (and mesh if provided)
            if (nargin >= 4)
                cOCPp.n_var     = varargin{1};
                cOCPp.n_contr   = varargin{2};
                cOCPp.xbc       = varargin{3};
                cOCPp.dyn_param = varargin{4};    
                if (nargin == 5)
                    cOCPp.mesh  = varargin{5};
                end
            end
        end

        
        % set methods
        function set.mesh(obj,mesh)
            % set method for the time mesh
            [n,m] = size(mesh);
            if ((n>0) && (m==1))
                obj.mesh = mesh';
            elseif ((n==1)&&(m>0))
                obj.mesh = mesh;
            else 
                error('Wrong dimensions in mesh input.');
            end
        end

        function set.xbc(obj,xbc)
            % set method for the boundary values
            [n,m] = size(xbc);
            if ((n>0) && (m==2))
                obj.xbc = xbc;
            elseif ((n==2)&&(m>0))
                obj.xbc = xbc';
            else 
                error('Wrong dimensions for boundary values (xbc) input.');
            end
        end        
        
        function set.n_var(obj,n_var)
            % set method for the number of state variables
            [n,m] = size(n_var);
            if ((n==1) && (m==1))
                if ((n_var > 0)&&(isinteger(n_var)))
                    obj.n_var = double(n_var);
                else 
                    error('n_var input is not a positive integer.');
                end
            else 
                error('Wrong dimensions in n_var input.');
            end
        end  
        
        function set.n_contr(obj,n_contr)
            % set method for the number of control variables
            [n,m] = size(n_contr);
            if ((n==1) && (m==1))
                if ((n_contr > 0)&&(isinteger(n_contr)))
                    obj.n_contr = double(n_contr);
                else 
                    error('n_var input is not a positive integer.');
                end
            else 
                error('Wrong dimensions in n_var input.');
            end
        end
        
        function set.dyn_param(obj,dyn_param)
            % set method for the classDynParam object
            if (isa(dyn_param,'classDynParam'))
                obj.dyn_param = dyn_param;
            else 
                error('Wrong type of dyn_param input.');
            end
        end       
        
        function set.n_intervals(obj,~)
            % the n_intervals property cannot be set but is a result of the
            % current mesh
            error('You cannot set n_intervals property.');
        end
        
        function set.n_timepoints(obj,~)
            % the n_timepoints property cannot be set but is a result of
            % the current mesh
            error('You cannot set n_timepoints property.');
        end
        
        % get routines
        function xbc = get.xbc(obj)
            % get the stored values for the boundary conditions
            [n,m] = size(obj.xbc);
            if ((n==obj.n_var)&&(m==2))
                xbc = obj.xbc;
            else
                error('Wrong dimensions in boundary values (xbc) output.');
            end
        end
        
        function n_intervals = get.n_intervals(obj)
            % get the number of intervals in the time mesh
            n_intervals = length(obj.mesh);
        end
        
        function n_timepoints = get.n_timepoints(obj)
            % get the number of time points in the mesh (n_intervals +1)
            n_timepoints = length(obj.mesh)+1;
        end        
        
        % other functions
        function setUniformMesh(obj,n_in)
            % function to set a uniform mesh with specified number of
            % intervals (without providing the actual mesh values)
            if ((isinteger(n_in))&&(n_in>1))
                obj.mesh = ones(1,n_in)*(1.0/(double(n_in)));
            else
                error('Input is no integer greater 1.');
            end
        end
          
    end
end