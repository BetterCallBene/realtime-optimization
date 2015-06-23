classdef TestInt < handle & TestEnv
    %TESTINT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dynForwEuler;
        dynRungeKutta;
        dynOde45;
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
            res = obj.dynRungeKutta.solver;
        end
        
        function res = get.solver3(obj)
            res = obj.dynOde45.solver;
        end
    end
    
    
    methods(Test)
        function testIntegratoren(testCase)
            n_intervals = 50;
            timepoint = 5;
            testCase.setupTest(n_intervals);
            [F1, J1] = testCase.solver1.odeTest(timepoint);
            [F2, J2] = testCase.solver2.odeTest(timepoint);
            [F3, J3] = testCase.solver3.odeTest(timepoint);
            

            anaDiff1 = J1;
            anaDiff2 = J2;
            
            func1 = @() testCase.solver1.odeTest(timepoint);
            numDiff1 = testCase.solver1.numDiff_nD(timepoint, func1);
            
            func2 = @() testCase.solver2.odeTest(timepoint);
            numDiff2 = testCase.solver2.numDiff_nD(timepoint, func2);
                
            testCase.assertLessThan(max(abs(anaDiff1 - anaDiff2)), 9e-2);
            
            %max(abs(anaDiff2 - numDiff2))
            %norm(anaDiff2 - numDiff2, 1)
        end
    end
    
    methods
        
        function dy = func(y)
            
        end
        
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
            
            
            %vec =  rand(17* (n_intervals+1), 1);
            %save('Data.mat', 'vec');
            vec = rand(17* (n_intervals+1), 1);
            
            solver1 = ForwEuler();
            solver2 = RungeKutta();
            solver3 = ode45M();
            
            obj.dynForwEuler = BasisQDyn(model, env, solver1);
            obj.dynForwEuler.vec = vec;
            
            obj.dynRungeKutta = BasisQDyn(model, env, solver2);
            obj.dynRungeKutta.vec = vec;
            
            obj.dynOde45 = BasisQDyn(model, env, solver3);
            obj.dynOde45.vec = vec;
        end
    end
    
end

