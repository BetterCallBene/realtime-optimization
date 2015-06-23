classdef ode45M < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods
        
        function M45 = ode45M()
            M45@Solver();
        end
        
        
        
        
        function y = integrate(obj, func, meshGrid, y0)
            [t, y] = ode45(func, meshGrid, y0);
            y = y(end, :)';
        end
        
        
    end
    
    methods(Test)
        
        function testOde(testCase)
            n_intervals = 50;
            
            testCase.setupTest(n_intervals);
            
            tic;
            
            for timepoint = 1:(n_intervals)
                
                [F, J, M, N] = testCase.odeTest(timepoint);
                                        
                func = @() testCase.odeTest(timepoint);
                numDiff = testCase.numDiff_nD(timepoint, func);
                testCase.assertLessThan(max(abs(J - numDiff)), testCase.tolRK);
            end
            
            toc
        end
        
    end
    
    
    
    
end

