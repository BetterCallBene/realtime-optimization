classdef Environment < handle
    % ENVIRONMENT Diese Klasse speichert alle Parameter, die den Quadrocopter nicht intern beschreiben, also Nebenbedingungen, Mesh, Gravitation, Wind (später)
    properties
        mesh;           % time mesh (1xn) vector
        xbc;            % (mx2) matrix with boundary values for states
        g = 9.81;		% Gravitation

        horizon;        % Horizon of the realtime approach
        wind;           % A function handle, which has s_t and t as input and returns the exact measured state
        n_intervals;     % number of intervals in mesh
    end
    
    properties(Dependent)
        n_timepoints;    % number of discretization points in mesh
    end
    
    methods
        % constructor
        function env = Environment(varargin)
            % constructor sets values for n_var, n_control, xbc and
            % dyn_param (and mesh if provided)
            if (nargin >= 1)
                env.xbc       = varargin{1};
                
                if (nargin == 2)
                    env.mesh  = varargin{2};
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
        
        % get routines
        function xbc = get.xbc(obj)
            % get the stored values for the boundary conditions
            [n,m] = size(obj.xbc);
            
            xbc = obj.xbc;
        end
        % 
%         function n_intervals = get.n_intervals(obj)
%             % get the number of intervals in the time mesh
%             n_intervals = length(obj.mesh);
%         end
        
       
        
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
                obj.n_intervals = length(obj.mesh);
            else
                error('Input is no integer greater 1.');
            end
        end
        
        function [n_intervals] = setUniformMesh1(obj,horizon, points_per_sec)
            % function to set a uniform mesh with specified number of
            % intervals (without providing the actual mesh values)
            if (horizon>1)
                obj.mesh = linspace(0, horizon, horizon * points_per_sec + 1);
                n_intervals = length(obj.mesh);
                obj.n_intervals = n_intervals;
            else
                error('Input is no integer greater 1.');
            end
%             if ((isinteger(n_in))&&(n_in>1))
%                 obj.mesh = ones(1,n_in)*(1.0/(double(n_in)));
%                 obj.n_intervals = length(obj.mesh);
%             else
%                 error('Input is no integer greater 1.');
%             end
            
        end
    end
end

