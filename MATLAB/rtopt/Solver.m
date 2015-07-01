classdef(Abstract) Solver < handle & TestEnv
    %SOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        opts;
        dyn_;
        timepoint; % 
        M0;
        N0;
        u0;
        %M0_vec;
        %N0_vec;
        M0_size;
        N0_size;
        
        n_intervalsInt;
        
        h;
        
        vec_sav;
        
        flagM0;
        flagN0;
        
        flagM0_size;
        flagN0_size
    end
    
    properties(Dependent)
        vec;
        dyn;
        JDot;
        JDot_u;
        JDot_x;
        contr;
    end
    
    
    methods(Abstract)
        integrate(obj, func, tspan, y0);
    end
    
    methods
        
        function set.dyn(obj, dyn)
            obj.dyn_ = dyn;
            obj.h =  obj.dyn.environment.mesh(1); 
        end
        
        function res =  get.dyn(obj)
            res = obj.dyn_;
        end
        
        function res = get.vec(obj)
            res = obj.dyn.backdoor_vec;
        end
        
        function set.vec(obj, vec)
            obj.dyn.backdoor_vec = vec;
        end
        
        function res = get.M0(obj)
             if obj.flagM0 == true
                 n_state = obj.dyn.robot.n_state;
                 obj.M0 = speye(n_state);
                 obj.flagM0 = false;
             end
             res = obj.M0;
         end
        
        function res = get.N0(obj)
            if obj.flagN0 == true
                [n_state, n_contr] = obj.getParams();
                obj.N0 = sparse(n_state, n_contr);
                obj.flagN0 = false;
            end
            res = obj.N0;
        end
        
        function res = get.M0_size(obj)
            
%             if obj.flagM0_size
%                 n_state = obj.getParams();
%                 obj.M0_size = n_state * n_state;
%                 obj.flagM0_size = false;
%             end
%             res = obj.M0_size;
            res =169; % Performance
        end
        function res = get.N0_size(obj)
%             if obj.flagN0_size
%                 [n_state, n_contr] = obj.getParams();
%                 obj.N0_size = n_state * n_contr;
%                 obj.flagN0_size = false;
%             end
%             res = obj.N0_size;
            res = 52; %Performance
        end
        
        function res = getJDot(obj)
            res = obj.dyn.dotD(obj.timepoint);
        end
        
        function res = get.JDot(obj)
            res = obj.getJDot();
        end
        function res = get.JDot_x(obj)
            n_state = obj.dyn.robot.n_state;
            res = obj.JDot(:, 1:n_state);
        end
        function res = get.JDot_u(obj)
            n_state = obj.dyn.robot.n_state;
            res = obj.JDot(:, n_state + 1:end);
        end
        
        function M = kM(obj, M0_)
           Jx = obj.JDot_x;
           M = Jx * M0_;
        end
        
        function N = kN(obj, N0_)
            Jx = obj.JDot_x;
            Ju = obj.JDot_u;
            N = Jx * N0_ +  Ju;
        end
        
    end
    
    methods %Help Functions
        function [n_state, n_contr, n_var] = getParams(obj)
            n_state =13; %= obj.dyn.robot.n_state;
            n_contr = 4;%= obj.dyn.robot.n_contr;
            n_var = 17;%n_var = obj.dyn.robot.n_var;
            %n_timepoints = obj.dyn.environment.n_timepoints;
        end
        function [F, M, N, J] = helperCreateMatrizen(obj, Y)
            
            %[n_state, n_contr] = obj.getParams();
            n_state = 13;
            n_contr = 4;
            
            M0_size_ = 169;%obj.M0_size; %Performance
            N0_size_ = 52;%obj.N0_size;
            
            F = Y(1:n_state, 1);
            
            M_vec = reshape(Y(n_state + 1: n_state + M0_size_, 1), [n_state, n_state]);
            M = sparse(M_vec);
            N_vec = reshape(Y(n_state + M0_size_ + 1: n_state + M0_size_ + N0_size_, 1), [n_state, n_contr]);
            N = sparse(N_vec);
            J = [M, N];
        end
        function y0 = helperCreateInitialConditions(obj, varargin)
            [n_state, n_contr, n_var] = obj.getParams();
            
            if (nargin == 2) % Bug in Matlab, es zaehlt obj mit
                obj.nextStep(varargin{1})
            end
                        
            %y0 = obj.helperCreateVektor(obj.dyn.state(:, obj.timepoint), obj.M0, obj.N0);
            state = obj.vec_sav((obj.timepoint - 1) * n_var +1:(obj.timepoint - 1) * n_var + n_state);
            y0 = obj.helperCreateVektor(state, obj.M0, obj.N0);
        end
        function res = helperCreateVektor(obj, F, M, N)
            %[n_state] = obj.getParams();
            n_state = 13;
            M0_size_ = obj.M0_size;
            N0_size_ = obj.N0_size;
            
            y = zeros(n_state + M0_size_ + N0_size_, 1);
            y(1:n_state, 1) = F;
            y(n_state + 1: n_state + M0_size_, 1) = reshape(M, [M0_size_, 1]);
            y(n_state + M0_size_ + 1: n_state + M0_size_ + N0_size_, 1) = reshape(N, [N0_size_, 1]);
            res = sparse(y);
        end
    end
    
    methods
        function s = Solver(n, varargin)
            
            if n == 1
                opts_ = varargin{1};
                s.opts = opts_{1};
            else
                s.opts = odeset('RelTol',1e-4,'AbsTol',1e-6);
            end
            
            s.flagM0_size = true;
            s.flagN0_size = true;
            s.flagM0 = true;
            s.flagN0 = true;
        end
        
        function val = get_contr(obj, timepoint)
           % (if not yet stored, extract them from vec)            
           
            [n_state, n_contr, n_var] = getParams(obj);

            val = obj.vec_sav((timepoint-1)*(n_var)+n_state + 1:...
                timepoint*(n_var));
            
        end
        
        function [old_intervals] = preToDo(obj)
            [n_state, n_contr, n_var] = getParams(obj);
            obj.vec_sav = obj.vec;  % Speichere alten vektor
            obj.vec = zeros(n_var, 1);
            old_intervals = obj.dyn.environment.n_intervals;
            obj.n_intervalsInt = old_intervals;
            obj.dyn.environment.n_intervals = 0;
        end
        
        function postToDo(obj, old_intervals)
            obj.vec = obj.vec_sav;
            obj.dyn.environment.n_intervals = old_intervals;
        end
        
        function [F, J, M, N] = ode(obj, timepoint, varargin)
            
            obj.nextStep(timepoint);
            
            if nargin <= 2
                y0 = obj.helperCreateInitialConditions();
                yp0 = [];
            elseif nargin == 4
                y0 = varargin{1};
                yp0 = varargin{2};
            else
                error('Wrong parameter count');
            end
            
            tspan = [(timepoint -1)*obj.h, timepoint*obj.h];
            meshGrid = [tspan(1), tspan(2)];%linspace(tspan(1), tspan(2));%, obj.n_intervalsInt); 
            
            old_timepoint = obj.timepoint;
            obj.u0 = obj.get_contr(old_timepoint);
            obj.timepoint = 1;
            
            y = obj.integrate(@obj.funcToIntegrate, meshGrid, y0, yp0);
            
            obj.timepoint = old_timepoint;
            [F, M, N, J] = obj.helperCreateMatrizen(y);
        end
        function [F, J, M, N] = odeTest(obj, timepoint, y0, yp0)
            [F, J, M, N] = obj.ode(timepoint, y0, yp0);
        end
        
        function dy = funcToIntegrate(obj, t, varargin)
            
            y = varargin{1};
            
            [n_state] = obj.getParams();
            
           
            [state, M0_, N0_] = obj.helperCreateMatrizen(y);            
            
            obj.vec = [state; obj.u0];
            
            dy = obj.helperCreateVektor(obj.dyn.dot(obj.timepoint), obj.kM(M0_), obj.kN(N0_));
        end
        
        function nextStep(obj, timepoint)
            obj.timepoint = timepoint;
        end
            

    end
    
    methods
        function setupTest(obj,n_intervals)
            % Quadrocopter soll 5 Meter hoch fliegen
            xbc = [         ... Variablenname L�nge   Name
                ... Anfangsbedingung
                0, 0, 0,    ...     r           3      Ortsvektor
                1, 0, 0, 0, ...     q           4      Quaternion (Einheitsquaternion)
                0, 0, 0,    ...     v           3      Translatorische Geschwindigkeit
                0, 0, 0;    ...     w           3      Winkelgeschwindigkeit
                ... Endbedingung
                0, 0, 5,    ...
                1, 0, 0, 0, ...
                0, 0, 0,    ...
                0, 0, 0     ...
                ];
            
            
            env = Environment();
            env.xbc = xbc;
            env.setUniformMesh(uint8(n_intervals));
                       
            model = Quadrocopter();
            
            obj.dyn = BasisQDyn(model, env, obj);
            obj.dyn.vec = rand(17* (n_intervals+1), 1);
            
        end
        
        function [vec_old, n, m, n_timepoints, dyn] = setup(obj,func)
            vec_old = obj.dyn.vec;
            n_timepoints = obj.dyn.environment.n_timepoints;
            n = obj.dyn.robot.n_var;
            dyn = obj.dyn;
            m = size(func());
            m=m(1);
        end
    end
    
    methods(Test)
    end
end

