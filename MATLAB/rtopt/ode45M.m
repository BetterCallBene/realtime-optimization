classdef ode45M < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    methods
        
        function M45 = ode45M(varargin)
            M45@Solver(nargin, varargin); %Bug in Matlab
        end
        
        
        function res = getJDot(obj)
            n_state = obj.dyn.robot.n_state;
            vec = obj.vec;
            state = vec(1:n_state);
            u0 = vec(n_state + 1 : end);
            res = obj.dyn.getJTilde(state, u0);
        end
        
        
        function y = integrate(obj, func, meshGrid, y0, yp0)
            opts_ = obj.opts;
            [t, y] = ode45(func, meshGrid, y0, opts_);
            y = y(end, :)';
        end
        
        function dy = funcToIntegrate(obj, t, varargin)
            
            y = varargin{1};
            u0 = obj.u0;
            
            [state, M0, N0] = obj.helperCreateMatrizen(y);
            
            obj.vec = [state; u0];
            FTilde = obj.dyn.FTilde(y, u0);
            
            kM = obj.kM(M0); 
            kN = obj.kN(N0);
            
            dy = obj.helperCreateVektor(FTilde, kM, kN);
        end
        
        function [y0, old_interval, old_timepoint] = setupTest(obj,n_intervals, timepoint)
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
            val = [zeros(3, 1); 1; 0 ;0 ;0; zeros(6, 1)];
            u = [10000; 10000; 10000+2000; 10000];
            vec = [val;u];
            obj.dyn.vec = [vec; vec; vec];
            
            
            obj.nextStep(timepoint);
            [old_interval] = obj.preToDo();
            
            y0 = obj.helperCreateInitialConditions();
            %yp0 = obj.helperCreateInitialConditionsDot();
            
            old_timepoint = obj.timepoint;
            obj.u0 = obj.get_contr(old_timepoint);
            obj.timepoint = 1;
            
        end
        
        
    end
    
    methods(Test)
        
        function testOde(testCase)
            n_intervals = 2;
            timepoint = 1;
            
            testCase.n_intervalsInt = n_intervals;
            
            [y0, old_interval, old_timepoint] = testCase.setupTest(n_intervals, timepoint);
            [n_state] = testCase.getParams();
            
            opts_ = odeset('RelTol',1e-6,'AbsTol',1e-7);%, 'Jacobian', @testCase.Jac);
            testCase.opts = opts_;
            
            tspan = [(timepoint -1)*testCase.h, timepoint*testCase.h];
            tic;
            [F] = testCase.odeTest(timepoint, y0, []);
            
            toc
            Q = norm(F(4:7));
            % Differenz zur 1
            testCase.assertLessThan(Q - 1,testCase.tolRK);
        end
        
    end
    
    
    
    
end

