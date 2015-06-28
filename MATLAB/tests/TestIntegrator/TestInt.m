classdef TestInt < handle & TestEnv
    %TESTINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dynForwEuler;
        dynOde15iM;
    end
    
    properties(Dependent)
        solver1;
        solver2;
        solver3;
    end
    
    methods
        function res = get.solver1(obj)
            res = obj.dynForwEuler.solver;
        end
        
        function res = get.solver2(obj)
            res = obj.dynOde15iM.solver;
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
            for i=1:n
                func_p = obj.plusEpsShift(i,timepoint,vec_old,func, n, solver);
                func_n = obj.minusEpsShift(i,timepoint,vec_old,func, n, solver);
                
                %Central difference
                numDiff(:,i) = (func_p - func_n)/2/obj.eps;
            end
        end
        
    end
    
    
    methods(Test)
        function testIntegratoren(testCase)
            n_intervals = 50;
            %timepoint = 5;
            
            testCase.setupTest(n_intervals);
%             [old_intervals] = testCase.solver1.preToDo();
%             [F1, J1] = testCase.solver1.ode(timepoint);
%             testCase.solver1.postToDo(old_intervals);
%             
%             
%             
%             [old_intervals] = testCase.solver1.preToDo();
%             func1 = @() testCase.solver1.ode(timepoint);
%             numDiff1 = testCase.numDiff_nD(timepoint, func1, testCase.solver1);
%             testCase.solver1.postToDo(old_intervals);
            
            for timepoint = 1:n_intervals
            
                [old_intervals] = testCase.solver2.preToDo();
                [F2, J2] = testCase.solver2.ode(timepoint);
                testCase.solver2.postToDo(old_intervals);
            
                [old_intervals] = testCase.solver2.preToDo();
                func2 = @() testCase.solver2.ode(timepoint);
                numDiff2 = testCase.numDiff_nD(timepoint, func2, testCase.solver2);
                testCase.solver2.postToDo(old_intervals);
                %figure
                timepoint
                max(max(abs(J2 - numDiff2)))
                
            end
            
            %func2 = @() testCase.solver2.odeTest(timepoint);
            %numDiff2 = testCase.solver2.numDiff_nD(timepoint, func2);
                
            %testCase.assertLessThan(max(abs(J1 - numDiff1)), 9e-2);
            
            %max(abs(anaDiff2 - numDiff2))
            %norm(anaDiff2 - numDiff2, 1)
        end
        
        %timepoint 5
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
            
            
            vec =  rand(17* (n_intervals+1), 1);
            %save('Data.mat', 'vec');
            
            %load('Data.mat', 'vec');
            
            
            solver1 = ForwEuler();
            solver2 = ode15iM2();
            
            
            obj.dynForwEuler = BasisQDyn(model, env, solver1);
            obj.dynForwEuler.vec = vec;
            
            obj.dynOde15iM = BasisQDyn(model, env, solver2);
            obj.dynOde15iM.vec = vec;
            
        end
    end
    
end

