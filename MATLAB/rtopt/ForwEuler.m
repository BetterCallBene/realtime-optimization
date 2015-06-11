classdef ForwEuler < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function FE = ForwEuler()
            FE@Solver(); %Bug von Matlab
        end
        
        function y = integrate(obj, func, tspan, y0)
            y = y0 + obj.h .* feval(func, y0);
        end
        
        function dy = funcToIntegrate(obj, y)
            
            [F, M, N] = obj.helperCreateMatrizen(y);
            
            M = obj.JDot_x * M;
            N = obj.JDot_x * N + obj.JDot_u;
            
            dy = helperCreateVektor(obj, obj.dyn.dot(obj.timepoint), M, N);
        end
        
        function [vec_old, n, m, n_timepoints, dyn] = setup(obj,func)
            vec_old = obj.vec;
            n_timepoints = obj.dyn.environment.n_timepoints;
            n = obj.dyn.robot.n_var;
            dyn = obj.dyn;
            m = size(func());
            m=m(1);
        end
        
    end
    
    methods(Test)
        function testIntegrate(testCase)
            n_intervals = 50;
            testCase.setupTest(n_intervals);
            tic;
            for timepoint = 1:(n_intervals)
            
                [F, J, M, N] = testCase.ode(timepoint);

                anaDiff = J;
            
                func = @() testCase.ode(timepoint);
                numDiff = testCase.numDiff_nD(timepoint, func);
                
                testCase.assertLessThan(max(abs(anaDiff - numDiff)), testCase.tol);
            end
            toc
        end
    end
    
    
    
    
end

