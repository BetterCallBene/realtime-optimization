classdef TestEnv < matlab.unittest.TestCase
    %TESTENV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        eps = 1e-2;
        tol = 1e-10;
    end
    
    methods
        
        function numDiff = numDiff_1D(obj,timepoint, func)
            % NUMDIFF1D This method calculates numerically the derivative numDiff
            % of func, when func only depends on obj.vec and has 1 dim output
            [vec_old, n,m] = obj.setup(func);
            numDiff = zeros(m,1);
            for i=1:n
                func_p = obj.plusEpsShift(i,timepoint,vec_old,func);
                func_n = obj.minusEpsShift(i,timepoint,vec_old,func);
                
                %Central difference
                numDiff(i) = (func_p - func_n)/2/obj.eps;
            end
        end
        
        function numDiff = numDiff_nD(obj,timepoint, func)
            % NUMDIFFDND This method calculates numerically the derivative
            % of func at timepoint t, when func only depends on obj.vec and has m dim output
            [vec_old, n,m] = obj.setup(func);
            numDiff = zeros(m,n);
            for i=1:n
                func_p = obj.plusEpsShift(i,timepoint,vec_old,func);
                func_n = obj.minusEpsShift(i,timepoint,vec_old,func);
                
                %Central difference
                numDiff(:,i) = (func_p - func_n)/2/obj.eps;
            end
        end
        
        function numDiff = numDiff_nD_AllT(obj, func)
            % NUMDIFF_NDALLT This method calculates numericalle the
            % derivative of func, but for all timepoints
            
            [vec_old, n, m, n_timestep] = obj.setup(func);
            numDiff = zeros(m * n_timestep, n * n_timestep);
            
            for timepoint = 1:n_timestep
                for i = 1:n
                    func_p = obj.plusEpsShift(i,timepoint,vec_old, func);
                    func_n = obj.minusEpsShift(i,timepoint,vec_old, func);
                    
                    %Central difference
                    numDiff((timepoint-1) * m : (timepoint*m) -1 , ((timepoint-1) * n) + i) ...
                        = (func_p - func_n)/2/obj.eps;
                end
            end
        end
        
        function numDiff = numDiff_nxnD(obj, timepoint, func)
            % NUMDIFFDNXND This method calculates numerically the
            % derivative of func, when func only depends on obj.vec and has
            % n times n dimensional output
            [vec_old, n, m] = obj.setup(func);
            numDiff = zeros(m,n,n);
            for i=1:n
                func_p = obj.plusEpsShift(i,timepoint, vec_old,func);
                func_n = obj.minusEpsShift(i,timepoint, vec_old,func);
                
                %Central difference
                numDiff(:,i,:) = (func_p - func_n )/2/obj.eps;
            end
        end
        
        
        %Some help functions
        function [vec_old, n, m, n_timepoints] = setup(obj, func)
            vec_old = obj.vec;
            n_timepoints = obj.environment.n_timepoints;
            n = obj.robot.n_var;
            m = size(func());
            m = m(1);
        end
        
        function func_p = plusEpsShift(obj,i,t,vec_old,func)
            vec_p = vec_old;
            vec_p((t-1)* obj.robot.n_var + i) = vec_p((t-1)* obj.robot.n_var + i) + obj.eps;
            obj.vec = vec_p;
            func_p = func();
            obj.vec = vec_old;
        end
        
        function func_n = minusEpsShift(obj,i,t,vec_old,func)
            vec_n = vec_old;
            vec_n((t-1)* obj.robot.n_var + i) = vec_n((t-1)* obj.robot.n_var + i) - obj.eps;
            obj.vec = vec_n;
            func_n = func();
            obj.vec = vec_old;
        end 
    end
end

