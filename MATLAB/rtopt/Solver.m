classdef(Abstract) Solver < handle & TestEnv
    %SOLVER Central providing the ODE solution
    
    
    properties
        opts;      % Speichert die ODE Options      
        dyn_;      % Handle auf die Klasse Dyn
        
        timepoint; % Aktueller Zeitpunkt
        M0;        % Initialisierung von M0 (Einheitsmatrix)
        N0;        % Initialisierung von M0 (Nullmatrix)
        
        
        M0_size;   % Groe�e der M0 - Matrix
        N0_size;   % Groe�e der N0 - Matrix
        
        
        h;         % Schrittlaenge
        u0;        % Aktuellen Controls
        vec_sav;   % Speicher des aktuellen vec ueber alle Zeitintervalle
        
        % Flags
        flagM0;
        flagN0;
        
        flagM0_size;
        flagN0_size
    end
    
    properties(Dependent)
        vec; % Steuerung des Zugriff auf den zentralen Vektor
        dyn; % Zugriff auf die Dynamik
        
        JDot;       % Aktuelle Jacobi der Dynamik zum Zeitpunkt timepoint
        JDot_x;     % Jacobimatrix der Dynamik ueber die States
        JDot_u;     % Jacobimatrix der Dynamik ueber die Controls
        
        contr;      % Controls
    end
    
    
    methods(Abstract)
        integrate(obj, func, tspan, y0, yp0); %Implementierung der Verfahren(Forward Euler, ode45..) in Subklassen
    end
    %Getter/Setter Methoden
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
            if obj.flagM0_size
                n_state = obj.getParams();
                obj.M0_size = n_state * n_state;
                obj.flagM0_size = false;
            end
            res = obj.M0_size;
            % Zur Performancesteigerung Optimierung hard coden
            %res =169; 
        end
        function res = get.N0_size(obj)
            if obj.flagN0_size
                [n_state, n_contr] = obj.getParams();
                obj.N0_size = n_state * n_contr;
                obj.flagN0_size = false;
            end
            res = obj.N0_size;
            % Zur Performancesteigerung Optimierung hard coden
            %res = 52; 
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
%     Help Functions
    methods 
        
        function [n_state, n_contr, n_var] = getParams(obj)
        % getParams Gibt folgende Parameter zurueck
        %   n_state, n_contr, n_var
        
            n_state = obj.dyn.robot.n_state;
            n_contr  =obj.dyn.robot.n_contr;
            n_var =  obj.dyn.robot.n_var;
            
%             Zur Performancesteigerung Optimierung hard coden
%             n_state =13; %= obj.dyn.robot.n_state;
%             n_contr = 4;%= obj.dyn.robot.n_contr;
%             n_var = 17;%n_var = obj.dyn.robot.n_var;
            
        end
        function [F, M, N, J] = helperCreateMatrizen(obj, Y)
        % helperCreateMatrizen Aus dem Initialisierungsvektor Y wird
        % folgendes generiert:
        %   F  enthaelt alle states und die Matrizen
        %   M  enthaelt die Jacobimatrix der states
        %   N  enthaelt die Jacobimatrix der controls
        %   J  gesamt Jacobimatrix
            
            [n_state, n_contr] = obj.getParams();
%             Zur Performancesteigerung Optimierung hard coden
%             n_state = 13;
%             n_contr = 4;
            
%             M0_size_ = 169;
%             N0_size_ = 52;
            
            M0_size_ = obj.M0_size; 
            N0_size_ = obj.N0_size;
            
            F = Y(1:n_state, 1);
            
            M_vec = reshape(Y(n_state + 1: n_state + M0_size_, 1), [n_state, n_state]);
            M = sparse(M_vec);
            N_vec = reshape(Y(n_state + M0_size_ + 1: n_state + M0_size_ + N0_size_, 1), [n_state, n_contr]);
            N = sparse(N_vec);
            J = [M, N];
        end
        function y0 = helperCreateInitialConditions(obj, varargin)
        % helperCreateInitialConditions Generiert aus den vec_sav und M0,
        % N0 den Initialisierungsvektor
       
            [n_state, n_contr, n_var] = obj.getParams();
            
            if (nargin == 2) % Bug in Matlab, es zaehlt obj mit
                obj.nextStep(varargin{1})
            end
            M0_ = obj.M0;
            N0_ = obj.N0;
            timepoint_ = obj.timepoint;
                        
            %y0 = obj.helperCreateVektor(obj.dyn.state(:, obj.timepoint), obj.M0, obj.N0);
            state = obj.vec_sav((timepoint_ - 1) * n_var +1:(timepoint_ - 1) * n_var + n_state);
            y0 = obj.helperCreateVektor(state, M0_, N0_);
        end
        function res = helperCreateVektor(obj, F, M, N)
        % helperCreateVektor Aus den Matrizen F, M und N wird ein langer
        % Vektor y erstellt
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
        
        function val = get_contr(obj, timepoint)
        % get_contr Extrahiert aus den Vektor obj.vec_sav die Controls zum
        % passend Zeitpunkt timepoint
           
            [n_state, n_contr, n_var] = getParams(obj);

            val = obj.vec_sav((timepoint-1)*(n_var)+n_state + 1:...
                timepoint*(n_var));
            
        end
        
        function [old_intervals] = preToDo(obj)
        % preToDo Initialisiert vor jedem Multishooting den Ode Solver
            [n_state, n_contr, n_var] = getParams(obj);
            obj.vec_sav = obj.vec;  % Speichere alten vektor
            obj.vec = zeros(n_var, 1);
            old_intervals = obj.dyn.environment.n_intervals;
            %obj.n_intervalsInt = old_intervals;
            obj.dyn.environment.n_intervals = 0;
        end
        
        function postToDo(obj, old_intervals)
        % postToDo Raeumt nach dem Multishooting auf
            obj.vec = obj.vec_sav;
            obj.dyn.environment.n_intervals = old_intervals;
        end
    end
    
    methods
        function s = Solver(n, varargin)
        % Solver Konstruktor initialisiert die Flags und setzt gebenfalls die Toleranz    
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
        
        function [F, J, M, N] = ode(obj, timepoint, varargin)
        % ode Fuehrt die Berechnung der Ode zum Zeitpunkt timepoint durch
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

            if((norm(y0(4:7)) - 1) > 1e-3)
                %warning('Quaternioen ist nicht normiert');
            end
            y = obj.integrate(@obj.funcToIntegrate, meshGrid, y0, yp0);
            
            if((norm(y(4:7)) - 1) > 1e-3)
                %warning('Quaternioen ist nicht normiert');
            end
            
            obj.timepoint = old_timepoint;
            [F, M, N, J] = obj.helperCreateMatrizen(y);
        end
        
        
        function dy = funcToIntegrate(obj, t, varargin)
        % funcToIntegrate Zu integrierende Funktion wird bei den meisten Subklassen ueberschrieben    
            y = varargin{1};
            u0_ = obj.u0;
            timepoint_ = obj.timepoint;
            
            [state, M0_, N0_] = obj.helperCreateMatrizen(y); 
            obj.vec = [state; u0_];
            
            F = obj.dyn.dot(timepoint_);
            M = obj.kM(M0_);
            N = obj.kN(N0_);
            dy = obj.helperCreateVektor(F, M, N);
        end
        
        function nextStep(obj, timepoint)
            obj.timepoint = timepoint;
        end
    end
    % Helper Funktion fuer die Tests
    methods
        function [F, J, M, N] = odeTest(obj, timepoint, y0, yp0)
            [F, J, M, N] = obj.ode(timepoint, y0, yp0);
        end
    end
    
    methods
        function setupTest(obj,n_intervals)
        %setupTest Initialisierung der Tests
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
        % setup Initialisierung der numerischen Ableitungen
            vec_old = obj.dyn.vec;
            n_timepoints = obj.dyn.environment.n_timepoints;
            n = obj.dyn.robot.n_var;
            dyn = obj.dyn;
            m = size(func());
            m=m(1);
        end
    end
end

