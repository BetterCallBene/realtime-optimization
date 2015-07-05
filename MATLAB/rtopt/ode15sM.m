classdef ode15sM < Solver
    %ODE15SM Integriert den MATLAB Solver ode15s in das Realtimeprojekt 
      
    
    properties
        y0;
    end
    % helper
    methods
        
%         function res = getJDot(obj)
%         % getJDot Verwendet die Funktion JTilde anstatt dotD in der Dynamik, um die
%         % Jacobimatrix zu bestimmen
%             n_state = obj.dyn.robot.n_state;
%             vec = obj.vec;
%             state = vec(1:n_state);
%             u0 = vec(n_state + 1 : end);
%             res = obj.dyn.getJTilde(state, u0);
%         end
        function yp = helperCreateInitialConditionsDot(obj)
        % helperCreateInitialConditionsDot Bestimmt die Initialisierung
        % Bedingung fuer die Ableitung zum Zeitpunkt timepoint
            [n_state, n_contr, n_var] = getParams(obj);
            
            old_timepoint = obj.timepoint;
            obj.timepoint = 1;
            vec = obj.vec_sav((old_timepoint - 1) * n_var +1:(old_timepoint - 1) * n_var + n_var);
            obj.dyn.vec = vec;
            
            xp = obj.dyn.dot(1);%obj.dyn.FTilde(vec(1:n_state, 1), vec(n_state+1:end));
            Mp = obj.kM(obj.M0);
            Np = obj.kN(obj.N0);
            
            obj.timepoint = old_timepoint;
            
            yp = obj.helperCreateVektor(xp, Mp, Np);
        end
        function dfdyX = JacX(obj, M0, N0, Jx, H)
        % JacX Bestimmt den state Anteil der Jacobimatrix
            [n_state, n_contr] = obj.getParams();
            
            Mx = zeros(n_state * n_state, n_state);
            Nx = zeros(n_state * n_contr, n_state);
            
            
            for i = 1:n_state
                M0spalte = full(M0(:, i));
                Hx = H(:, 1:n_state, 1:n_state);
                Mx((i-1)*n_state+1:i*n_state, :) = tprod(Hx, [1 -1 2], M0spalte, [-1]);
                if i <= n_contr
                    N0spalte = full(N0(:, i));
                    
                    Nx((i-1)*n_state+1:i*n_state, :) = tprod(Hx, [1 -1 2], N0spalte, [-1])  ...
                        + tprod(1, [-1],  H(:, i + n_state, 1:n_state), [1 -1, 2]);
                end
            end
            
            dfdyX = [sparse(Jx);
                sparse(Mx);
                sparse(Nx);
            ];
            
        end
        
        function dfdyM = JacM(obj, Jx)
        % JacM Bestimmt den M Anteil der Jacobimatrix
            Mdiag = Jx;
            
            dfdyM = [
                    sparse(13, 169);
                    sparse(blkdiag(Mdiag, Mdiag, Mdiag, Mdiag, Mdiag,...
                Mdiag, Mdiag, Mdiag, Mdiag, Mdiag,....
                Mdiag, Mdiag, Mdiag ...
            ));
                sparse(4*13, 169);
            ];
            
        end
        
        function dfdyM = JacN(obj, Jx)
        % JacN Bestimmt den M Anteil der Jacobimatrix   
            Ndiag = Jx;
            
            dfdyM = [
                    sparse(13, 52);
                    sparse(169, 52);
                    sparse(blkdiag(Ndiag, Ndiag, Ndiag, Ndiag));
            ];
            
        end
    end
    methods
        
        function M15s = ode15sM(varargin)
        % ode15sM Ruft den Konstruktor Basisklasse auf und setzt den
        % Funktion handler fuer die Jacobimatrix des ode15s Solvers
            M15s@Solver(nargin, varargin); %Bug in Matlab
            M15s.opts.Jacobian = @M15s.Jac; 
        end
        
        
        
        
        function y = integrate(obj, func, meshGrid, y0, yp0)
        % integrate Initialisert den Solver ode15s und fuehrt in aus    
            opts_ = obj.opts;
            opts_.InitialSlope = yp0;
            
            [t, y] = ode15s(func, meshGrid, y0, opts_);
            y = y(end, :)';
        end
        
        function [F, J, M, N] = ode(obj, timepoint, varargin)
        % ode Erweitert die Funktion ode in Basisklasse mit den
        % Initialisierungsbedingungen der Ableitung zum Zeitpunkt timepoint
            obj.nextStep(timepoint);
            
            y0_ = obj.helperCreateInitialConditions();
            yp0 = obj.helperCreateInitialConditionsDot();
            [F, J, M, N] = ode@Solver(obj, timepoint, y0_, yp0);
        end
        
        function dy = funcToIntegrate(obj, t, varargin)
        % funcToIntegrate Verwendet die Funktion FTilde aus Dynamik anstatt
        % ode
            y = varargin{1};
            u0 = obj.u0;
            
            [state, M0, N0] = obj.helperCreateMatrizen(y);
            
            obj.vec = [state; u0];
            F = obj.dyn.dot(1);
            
            kM = obj.kM(M0); 
            kN = obj.kN(N0);
            
            dy = obj.helperCreateVektor(F, kM, kN);
        end
        
        function dfdy = Jac(obj, t, y)
        % Jac Bestimmt die Jacobimatrix fuer das ode15s Verfahren    
            u0 = obj.u0;
            [x0, M0, N0] = obj.helperCreateMatrizen(y);
            
            obj.vec = [x0; u0];
            
            %JTilde = obj.dyn.getJTilde(x0, u0);
            H = obj.dyn.dotDD(1);%getHTilde(x0, u0);
            Jx = obj.JDot_x; %JTilde(1:n_state, 1:n_state);
            
            dfdy = [obj.JacX(M0, N0, Jx, H), ...
                    obj.JacM(Jx), ...
                    obj.JacN(Jx) ...
                    ];
            
        end
        
    end
    %Helper fuer Testfunktionen
    methods
        
        function [vec_old, n, m, timepoints, dyn]  = setup(obj,func)
            y0_ = obj.y0;
            
            old_timepoint = obj.timepoint;
            %obj.u0 = obj.get_contr(old_timepoint);
            
            obj.nextStep(old_timepoint);
            
            timepoints = 0;
            dyn = obj.dyn;
            
            %n = 17;
            m = size(func([], y0_), 1);
            n = m;
            
            vec_old = y0_;
            
        end
        
        function  setupTest(obj,n_intervals, timepoint)
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
            obj.u0 = [10000; 10000; 10000; 10000];
            obj.y0 = rand(234, 1);
            
        end
        
        function func_p = plusEpsShift(obj,i,t,vec_old,func, n_var, dyn)
            %vec_old = dyn.vec;
            vec_p = vec_old;
            t = 1;
            vec_p((t-1)* n_var + i) = vec_p((t-1)* n_var + i) + obj.eps;
            %dyn.backdoor_vec = vec_p;
            func_p = func([], vec_p);
            
        end
        
        function func_n = minusEpsShift(obj,i,t,vec_old,func, n_var, dyn)
            %vec_old = dyn.vec;
            vec_n = vec_old;
            t = 1;
            vec_n((t-1)* n_var + i) = vec_n((t-1)* n_var + i) - obj.eps;
            func_n = func([], vec_n);
        end
    end
    
    methods(Test)
        
        function testJac(testCase)
        % testJac Test, ob die Jac Matrix zu der Funktion funcToIntegrate passt
            n_intervals = 4;
            timepoint = 1;
            testCase.timepoint = timepoint;
            testCase.setupTest(n_intervals, timepoint);
            [old_interval] = testCase.preToDo();
            [anaJ] =testCase.Jac([], testCase.y0);
            testCase.postToDo(old_interval);
            
            [old_interval] = testCase.preToDo();
            numDiffJ = testCase.numDiff_nD(timepoint, @testCase.funcToIntegrate);
            testCase.postToDo(old_interval);
            testCase.assertLessThan(max(abs(anaJ - numDiffJ)),testCase.tol);
        end
        
    end
    
    
    
    
end

