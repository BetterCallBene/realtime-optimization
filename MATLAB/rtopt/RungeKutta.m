classdef RungeKutta < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        u0;
    end
    
    methods
        
        function RK = RungeKutta()
            RK@Solver();
        end
        
        function M = kM(obj, y)
            n_state  = obj.dyn.robot.n_state;
            M0 = reshape(y, [n_state, n_state]);
            M = obj.JDot_x * M0;
        end
        
        function dy = funcToIntegrate(obj, t, y)
            [n_state] = obj.getParams();
            
            timepoint = 1;
            
            obj.vec = [y(1:n_state); obj.u0];
            dy = obj.helperCreateVektor(obj.dyn.dot(timepoint), obj.kM(y(n_state+1:n_state + obj.M0_size)), obj.JDot_u);
        end
        
        function y = RungeKuttaSolver(obj, func, meshGrid, y0)
            h =  obj.h;
            y = y0;
            for i = 1:1
                k1 = func([], y);
                k2 = func([], y + h/2.*k1);
                k3 = func([], y + h/2.*k2);
                k4 = func([], y + h .* k3);
                y =  y +h/6.*(k1 + 2.* k2 + 2.* k3 + k4);
            end
        end
        
        function [old_vec, old_timepoint, old_intervals] = preToDo(obj)
            timepoint = obj.timepoint; 
            obj.u0 = obj.dyn.contr(:, timepoint);
            
            old_vec = obj.vec;  % Speichere alten vektor
            old_timepoint = timepoint;
            old_intervals = obj.dyn.environment.n_intervals;
            
            obj.timepoint = 1;
            obj.dyn.environment.n_intervals = 0;
        end
        
        function postToDo(obj, old_vec, old_timepoint, old_intervals)
            obj.dyn.environment.n_intervals = old_intervals;
            obj.timepoint = old_timepoint;
            obj.vec = old_vec;
        end
        
        function Y = integrate(obj, func, meshGrid, y0)
            
            [old_vec, old_timepoint, old_intervals] = obj.preToDo();
            [Y] = obj.RungeKuttaSolver(func, meshGrid, y0);%ode45(@obj.kFDot_Test, meshGrid, y0);
            obj.postToDo(old_vec, old_timepoint, old_intervals);
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
        
        function testOde(testCase)
            n_intervals = 50;
            
            testCase.setupTest(n_intervals);
            
            tic;            
            for timepoint = 1:(n_intervals)
                [F, J, M, N] = testCase.ode(timepoint);
                                        
                func = @() testCase.ode(timepoint);
                numDiff = testCase.numDiff_nD(timepoint, func);
                testCase.assertLessThan(max(abs(J - numDiff)), testCase.tolRK);
            end
            toc
        end
        
    end
    
    
    
    
end

