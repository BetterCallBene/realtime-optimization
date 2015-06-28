classdef TestInt < handle & TestEnv
    %TESTINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dynOde45M;
        dynOde15iM;
        dynOde15sM;
    end
    
    properties(Dependent)
        solvers;
        solver1;
        solver2;
        solver3;
    end
    
    methods
        function res = get.solver1(obj)
            res = obj.dynOde45M.solver;
        end
        
        function res = get.solver2(obj)
            res = obj.dynOde15iM.solver;
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
        
        function [error, norm_error] = NumAnaTest(obj, n_intervals, solver)
            error = zeros(n_intervals, 1);
            norm_error = zeros(n_intervals, 1);
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
                norm_error(timepoint) = norm(F(4:7)) - 1;
            end
        end
        
    end
    
    
    methods(Test)
        function SpeedTestSolvers(testCase)
            n_intervals = 20;
           
            testCase.setupTest(n_intervals);
            
            solver = testCase.solver3;
            
            
            %n_intervals = 50;
            [old_intervals] = solver.preToDo();
            tic
            parfor timepoint = 1:n_intervals
                [F3, J3] = solver.ode(timepoint);
            end
            toc
            solver.postToDo(old_intervals);
            
        end
        function testIntegratoren(testCase)
            n_intervals = 20;
            
            %testCase.setupTest(n_intervals);
            testCase.setupTest2(n_intervals);
            % Fuer qualitative Aussagen mehrmals solver ausführen
            [spd] = testCase.speedTest(testCase.solver1, n_intervals);
            [spd] = testCase.speedTest(testCase.solver1, n_intervals);
            [spd] = testCase.speedTest(testCase.solver1, n_intervals);
            [spd] = testCase.speedTest(testCase.solver1, n_intervals);
            [spd] = testCase.speedTest(testCase.solver1, n_intervals);
            disp('SpeedTest: Solver1');
            [spd] = testCase.speedTest(testCase.solver1, n_intervals);
            disp(spd)
            disp('SpeedTest: Solver3');
            [spd] = testCase.speedTest(testCase.solver3, n_intervals);
            disp(spd);
            pause
            disp('Numerikvergleich von Solver1');
            [error, norm_error] =  testCase.NumAnaTest(n_intervals, testCase.solver1);
            disp('Error')
            disp(error)
            disp('Norm error')
            pause
            disp(norm_error)
            pause
            disp('Numerikvergleich von Solver3');
            [error, norm_error] = testCase.NumAnaTest(n_intervals, testCase.solver3);
            disp('Error')
            disp(error)
            pause
            disp('Norm error')
            disp(norm_error)
            
        end
        
        %timepoint 5
    end
    
    methods
        
        function setupTest(obj,n_intervals)
            
            % Quadrocopter soll 5 Meter hoch fliegen
            xbc = [         ... Variablenname Lï¿½nge   Name
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
            
            opts_ = odeset('RelTol',1e-4,'AbsTol',1e-6);
            
            solver1 = ode45M(opts_);%ForwEuler();
            solver2 = ode15iM2(opts_);
            solver3 = ode15sM(opts_);
            
            
            obj.dynOde45M = BasisQDyn(model, env, solver1);
            obj.dynOde45M.vec = vec;
            
            obj.dynOde15iM = BasisQDyn(model, env, solver2);
            obj.dynOde15iM.vec = vec;
            
            obj.dynOde15sM = BasisQDyn(model, env, solver3);
            obj.dynOde15sM.vec = vec;
            
        end
        function setupTest2(obj,n_intervals)
            
            % Quadrocopter soll 5 Meter hoch fliegen
            xbc = [         ... Variablenname Lï¿½nge   Name
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
            
            
            val = [zeros(3, 1); 0; 0 ;0 ;1; zeros(6, 1)];
            u = [10000, 10000, 10000, 10000+3000];
            vec = [val; u'];
            
            vec = repmat(vec, (n_intervals + 1), 1);
            %vec =rand(17 * (n_intervals + 1), 1); 
            
            opts_ = odeset('RelTol',1e-5,'AbsTol',1e-7);
            
            solver1 = ode45M(opts_);%ForwEuler();
            solver2 = ode15iM2(opts_);
            solver3 = ode15sM(opts_);
            
            
            obj.dynOde45M = BasisQDyn(model, env, solver1);
            obj.dynOde45M.vec = vec;
            
            obj.dynOde15iM = BasisQDyn(model, env, solver2);
            obj.dynOde15iM.vec = vec;
            
            obj.dynOde15sM = BasisQDyn(model, env, solver3);
            obj.dynOde15sM.vec = vec;
            
        end
    end
    
end

