classdef Constraints < GenConstraints & TestEnv
% classConstraints providing discretized ODE constraint using forward euler

    properties
        
    end
    
    methods
        %constructor
        function cC = Constraints(varargin)
            dode = [];
            dyn = [];
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
            elseif (nargin == 1)
                if (isa(varargin{1},'ForwEuler'))
                    dode = varargin{1};
                    dyn = cGC.dode.dyn; % Initializierung von Dynamik
                else
                    error('wrong class type for discretized ode');
                end
            else
                error('wrong number of inputs');     
            end
            cC@GenConstraints(dode, dyn);
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
        %BB: Nebenbedingung: (Norm(q))^2 = 1 hinzugefügt
        function eq_con = get_eq_con(obj)
            % the equality constraint of the ocp
            % combine the discretized ode with the boundary conditions
            
            xbc         = obj.dyn.environment.xbc;
            state_mat   = obj.dyn.state;
            
            eq_con      = [obj.dode.h(); 
                                state_mat(:,1) - xbc(:,1); ...
                                state_mat(:,end) - xbc(:,end); ...
                                obj.EqCon ...
                          ];
        end
        
        function eq_conD = get_eq_conD(obj)
            % the Jacobian of the equality contraints of the ocp
            [q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);

            %q           = state_mat(4:7, :);
            
            srow        = 1:2*n_state;
            scol        = [1:n_state,...
                n_int*(n_state+n_contr)+1:n_int*(n_state+n_contr)+n_state];
            sval        = ones(1,2*n_state);

            eq_conD     = [obj.dode.hD(); sparse(srow,scol,sval,...
                                2*n_state,(n_int+1)*(n_state+n_contr)); ...
                                obj.EqConD ...
                          ]';

        end
        
        function eq_conDD = get_eq_conDD(varargin)
            % the Hessian of the equality constraints of the ocp.
            % Returns a cell array of Hessians if only the object 
            % is provided. If the object and the Langrange multipliers are
            % provided, the function returns the Hessian of the Lagrangian
            if (nargin == 1)
                obj         = varargin{1};
                [q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
            
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
    
    methods(Test)
        function test_get_eq_conD(obj)
            n_int_ = uint16(50);
            % Quadrocopter soll 5 Meter hoch fliegen
            xbc = [         ... Variablenname Lï¿½nge   Name
                ... Anfangsbedingung
                0, 0, 0.5,  ...     r           3      Ortsvektor
                1, 0, 0, 0, ...     q           4      Quaternion (Einheitsquaternion)
                0, 0, 0,    ...     v           3      Translatorische Geschwindigkeit
                0, 0, 0;    ...     w           3      Winkelgeschwindigkeit
                ... Endbedingung
                0, 0, 0,    ...
                1, 0, 0, 0, ...
                0, 0, 0,    ...
                0, 0, 0     ...
                ];    

            env = Environment();
            env.xbc = xbc;
            env.setUniformMesh(n_int_);
            
            robot = Quadrocopter();
            n_state_ = robot.n_state;
            n_contr_ = robot.n_contr;
            
            n_var = n_state_ + n_contr_;
            
            dyn_ = BasisQDyn(robot, env);
            dyn_.vec = rand(n_var * (n_int_ + 1), 1);
            
            obj.dode = ForwEuler(dyn_);
            obj.dyn = obj.dode.dyn;
            
            val0 = obj.get_eq_con();
            val1 = obj.get_eq_conD();
        end
    end
end