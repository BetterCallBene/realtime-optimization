classdef ForwEuler < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        function FE = ForwEuler()
            FE@Solver(); %Bug von Matlab
        end
        
        function y = integrate(obj, func, meshGrid, y0)
            y = y0 + obj.h .* feval(func, meshGrid, y0);
        end
        
    end
    
    methods(Test)
        function testIntegrate(testCase)
            n_intervals = 50;
            testCase.setupTest(n_intervals);
            tic;
            for timepoint = 1:(n_intervals)
                
                [F, J, M, N] = testCase.odeTest(timepoint);

                anaDiff = J;
            
                func = @() testCase.odeTest(timepoint);
                numDiff = testCase.numDiff_nD(timepoint, func);
                
                testCase.assertLessThan(max(abs(anaDiff - numDiff)), testCase.tol);
            end
            toc
        end
    end
    
    
    
    
end

