classdef TestInt < handle & TestEnv
    %TESTINT Testklasse bzw. Testumgebung der Integrator
    %   In dieser Klasse werden die Solver ode45, ode15s und Einschritt
    %   ForwEuler auf Geschwindigkeit sowie Genauigkeit getestet. Zudem
    %   werden die Sensitivitaeten getestet
    
    properties
        dynOde45M; %Dynamikklasse, die mit ode45 Verfahren initialisiert wird
        dynForwEuler; % " ForwEuler "
        dynOde15sM;  % " ode15s "
    end
    
    properties(Dependent)
        solver1; % Shortcout fuer ForwEuler
        solver2; % " ode45
        solver3; % " ode15s
    end
    
    methods
        
        function res = get.solver1(obj)
            res = obj.dynForwEuler.solver;
        end
        function res = get.solver2(obj)
            res = obj.dynOde45M.solver;
        end
        
        
        function res = get.solver3(obj)
            res = obj.dynOde15sM.solver;
        end
        
        
        %Some help functions (typically overwritten in subclasses)
        function [vec_old, n, m, n_timepoints] = setup(obj, func, solver)
            vec_old = solver.vec_sav;
            n_timepoints = solver.dyn.environment.n_timepoints;
            n = solver.dyn.robot.n_var;
            m = size(func());
            m = m(1);
        end
        
        %Overwrite function in TestEnv, because we don't want normed
        %quaternios here.
        function func_p = plusEpsShift(obj,i,t,vec_old,func, n_var, solver)
            %vec_old = dyn.vec;
            vec_p = vec_old;
            vec_p((t-1)* n_var + i) = vec_p((t-1)* n_var + i) + obj.eps;
            solver.vec_sav = vec_p;
            func_p = func();
            solver.vec_sav = vec_old;
        end
        
        %Overwrite function in TestEnv, because we don't want normed
        %quaternios here.
        function func_n = minusEpsShift(obj,i,t,vec_old,func, n_var, solver)
            %vec_old = dyn.vec;
            vec_n = vec_old;
            vec_n((t-1)* n_var + i) = vec_n((t-1)* n_var + i) - obj.eps;
            solver.vec_sav = vec_n;
            func_n = func();
            solver.vec_sav = vec_old;
        end
        
        function numDiff = numDiff_nD(obj, timepoint, func, solver)
            % NUMDIFFDND This method calculates numerically the derivative
            % of func at timepoint t, when func only depends on obj.vec and has m dim output
            [vec_old, n, m, timepoints] = obj.setup(func,solver);
            numDiff = zeros(m,n);
            parfor i=1:n
                func_p = obj.plusEpsShift(i,timepoint,vec_old,func, n, solver);
                func_n = obj.minusEpsShift(i,timepoint,vec_old,func, n, solver);
                
                %Central difference
                numDiff(:,i) = (func_p - func_n)/2/obj.eps;
            end
        end
        
        function [spd] = speedTest(obj, solver, n_intervals)
            
            [old_intervals] = solver.preToDo();
            tic;
            parfor timepoint = 1:n_intervals
                [F, J] = solver.ode(timepoint);
            end
            spd = toc;
            solver.postToDo(old_intervals);
        end
        
        function [error, norm_error, rel_error] = NumAnaTest(obj, n_intervals, solver)
            %NumAnaTest Test, die Sensitivitaeten auf Korrektheit
            error = zeros(n_intervals, 1);
            norm_error = zeros(n_intervals, 1);
            rel_error = zeros(n_intervals, 1);
            parfor timepoint = 1:n_intervals
                [old_intervals] = solver.preToDo();
                [F, J] = solver.ode(timepoint);
                solver.postToDo(old_intervals);

                [old_intervals] = solver.preToDo();
                func = @() solver.ode(timepoint);
                numDiff = obj.numDiff_nD(timepoint, func, solver);
                solver.postToDo(old_intervals);
                
                %error(max(max(abs(J - numDiff))));
                error(timepoint) = norm(J - numDiff, Inf);
                rel_error(timepoint) = error(timepoint) / norm(J, Inf);
                norm_error(timepoint) = abs(norm(F(4:7)) - 1);
            end
        end
        
        function stressTestIntegrator(testCase, n_intervals)
            % stressTestIntegrator Stresstest fuer die Integratoren
            
            solver1 =testCase.solver1;
            opts1 = solver1.opts;
            
            solver2 =testCase.solver2;
            opts2 = solver2.opts;
            
            solver3 =testCase.solver3;
            opts3 = solver3.opts;
            
            % Fuer qualitative Aussagen mehrmals solver ausf�hren
            [spd] = testCase.speedTest(solver1, n_intervals);
            [spd] = testCase.speedTest(solver1, n_intervals);
            [spd] = testCase.speedTest(solver1, n_intervals);
            [spd] = testCase.speedTest(solver1, n_intervals);
            [spd] = testCase.speedTest(solver1, n_intervals);
            disp('SpeedTest: ForEuler');
            [spd] = testCase.speedTest(solver1, n_intervals);
            disp(spd)
            disp('SpeedTest: ode45M');
            [spd] = testCase.speedTest(solver2, n_intervals);
            disp(spd);
            disp('SpeedTest: ode15sM');
            [spd] = testCase.speedTest(solver3, n_intervals);
            disp(spd);
            pause(1)
            disp('Numerikvergleich von ForwEuler');
            [error, norm_error, rel_error1] =  testCase.NumAnaTest(n_intervals, solver1);
            testCase.verifyLessThan(max(error), opts1.RelTol * 150);
            testCase.verifyLessThan(max(norm_error), 1e-3);
            
            disp('Numerikvergleich von ode45');
            [error, norm_error, rel_error2] =  testCase.NumAnaTest(n_intervals, solver2);
            testCase.verifyLessThan(max(error), opts1.RelTol * 150);
            testCase.verifyLessThan(max(norm_error), 1e-3);
            
            disp('Numerikvergleich von ode15sM');
            [error, norm_error, rel_error3] = testCase.NumAnaTest(n_intervals, solver3);
            testCase.verifyLessThan(max(error), opts3.RelTol * 150);
            testCase.verifyLessThan(max(norm_error), opts3.RelTol);
        end
        
        
         
    end
    
    
    methods(Test)
        function testSetup1(testCase)
            %testSetup1 Test der Integratoren mit Randomwerten
            n_intervals = 50;
            
            %testCase.setupTest(n_intervals);
            testCase.setupTest(n_intervals);
            disp('TestInt: TestSetup 1');
            
            testCase.stressTestIntegrator(n_intervals);
            
        end
        
        function testSetup2(testCase)
            %testSetup2 Test der Integratoren mit dem Punkt [1 0 0 0]
            n_intervals = 50;
            
            %testCase.setupTest(n_intervals);
            
            testCase.setupTest2(n_intervals);
            disp('TestInt: TestSetup 2');
            testCase.stressTestIntegrator(n_intervals);
        end
        
        function testExactness(testCase)
            %testExactness Vergleichstest der einzelnen Solvers in Bezug zu
            %ihrer Genauigkeit
            n_intervals = 50;
            testCase.setupTest(n_intervals)
            solver1 = testCase.solver1;
            solver2 = testCase.solver2;
            solver3 = testCase.solver3;
            
            solver1.preToDo();
            solver2.preToDo();
            [old_intervals] = solver3.preToDo();
            
            res = zeros(n_intervals, 3);
            tic;
            for timepoint = 1:n_intervals
                [F1, J1] = solver1.ode(timepoint);
                [F2, J2] = solver2.ode(timepoint);
                [F3, J3] = solver3.ode(timepoint);
                res(timepoint, 1) = max(abs(F1 - F2));
                res(timepoint, 2) = max(abs(F2 - F3));
                res(timepoint, 3) = max(abs(F3 - F1));
            end
            
            mittelF1F2 = sum(res(:, 1))/n_intervals;
            mittelF2F3 = sum(res(:, 2))/n_intervals;
            mittelF3F1 = sum(res(:, 3))/n_intervals;
            
            maxF1F2 = max(res(:, 1));
            maxF2F3 = max(res(:, 2));
            maxF3F1 = max(res(:, 3));
            
            disp(['Vergleich im Mittel zwischen Euler und Ode45: ', num2str(mittelF1F2)]);
            disp(['Vergleich im Mittel zwischen ode15s und Ode45: ', num2str(mittelF2F3)]);
            disp(['Vergleich im Mittel zwischen Euler und Ode45: ', num2str(mittelF3F1)]);
            pause(5)
            disp('');
            disp(['Maximalerunterschied zwischen Euler und Ode45: ', num2str(maxF1F2)]);
            disp(['Maximalerunterschied zwischen ode15s und Ode45: ', num2str(maxF2F3)]);
            disp(['Maximalerunterschied zwischen Euler und Ode45: ', num2str(maxF3F1)]);
            spd = toc;
            solver1.postToDo(old_intervals);
            solver2.postToDo(old_intervals);
            solver3.postToDo(old_intervals);
            
        end
        
        %timepoint 5
    end
    
    methods % Verschiedene Konfigurationen
        
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
            
            
            %val = [zeros(3, 1); 1; 0 ;0 ;0; zeros(6, 1)];
            %u = [10000, 10000, 10000, 10000+400];
            %vec = [val; u'];
            %vec = [vec; vec; vec];
            vec =rand(17 * (n_intervals + 1), 1); 
            
            opts_ = odeset('RelTol',1e-7,'AbsTol',1e-8);
            
            solver1 = ForwEuler();
            solver2 = ode45M(opts_);
            solver3 = ode15sM(opts_);
            
            obj.dynForwEuler = BasisQDyn(model, env, solver1);
            obj.dynForwEuler.vec = vec;
            
            obj.dynOde45M = BasisQDyn(model, env, solver2);
            obj.dynOde45M.vec = vec;
          
            obj.dynOde15sM = BasisQDyn(model, env, solver3);
            obj.dynOde15sM.vec = vec;
            
        end
        function setupTest2(obj,n_intervals)
            
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
            
            
            val = [zeros(3, 1); 1; 0 ;0 ;0; zeros(6, 1)];
            u = [5, 5, 5, 5];
            vec = [val; u'];
            
            vec = repmat(vec, (n_intervals + 1), 1);
            %vec =rand(17 * (n_intervals + 1), 1); 
            
            opts_ = odeset('RelTol',1e-4,'AbsTol',1e-5);
            
            solver1 = ForwEuler();
            solver2 = ode45M(opts_);
            solver3 = ode15sM(opts_);
            
            
            
            obj.dynForwEuler = BasisQDyn(model, env, solver1);
            obj.dynForwEuler.vec = vec;
            
            obj.dynOde45M = BasisQDyn(model, env, solver2);
            obj.dynOde45M.vec = vec;
            
            obj.dynOde15sM = BasisQDyn(model, env, solver3);
            obj.dynOde15sM.vec = vec;
            
        end
    end
    
end

