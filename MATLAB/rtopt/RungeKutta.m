classdef RungeKutta < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods
        
        function RK = RungeKutta()
            RK@Solver();
        end
        
        
        
        
        function y = integrate(obj, func, meshGrid, y0)
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

