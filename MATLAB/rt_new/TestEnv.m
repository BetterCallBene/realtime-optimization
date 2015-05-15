classdef TestEnv < matlab.unittest.TestCase
    %TESTENV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        eps = 1e-2;
    end
    
    methods
        
        %TODO: debug
        function numDiff = numDiff1D(obj,func)
            % NUMDIFF1D This method calculates numerically the derivative numDiff
            % of func, when func only depends on obj.vec and has 1 dim output
            [vec_old, n,m] = obj.setup(g_obj);
            numDiff = zeros(n,1);
            for i=1:n
                func_p = obj.plusEpsShift(i,vec_old,func);
                func_n = obj.minusEpsShift(i,vec_old,func);
                
                %Central difference
                numDiff(i) = (func_p - func_n)/2/obj.eps;
            end
        end
        
        %TODO: debug
        function numDiff = numDiffDnD(obj, func)
            % NUMDIFFDND This method calculates numerically the derivative
            % of func, when func only depends on obj.vec and has m dim output
            [vec_old, n,m] = obj.setup(g_obj);
            numDiff = zeros(m,n);
            for i=1:n
                func_p = obj.plusEpsShift(i,vec_old,func);
                func_n = obj.minusEpsShift(i,vec_old,func);
                
                %Central difference
                numDiff(:,i) = (func_p - func_n)/2/obj.eps;
            end
        end
        
        %TODO: debug
        function numDiff = numdiffDnxnD(obj,func)
            % NUMDIFFDNXND This method calculates numerically the
            % derivative of func, when func only depends on obj.vec and has
            % n times n dimensional output
            [vec_old, n] = obj.setup(g_obj);
            
            numDiff = zeros(m,n,n);
            for i=1:n
                func_p = obj.plusEpsShift(i,vec_old,func);
                func_n = obj.minusEpsShift(i,vec_old,func);
                
                %Central difference
                numDiff(:,i,:) = (func_p - func_n )/2/obj.eps; 
            end
        end
        
        function [vec_old, n,m] = setup(obj,g_obj,func)
            vec_old = g_obj.vec;
            n = length(vec_old);
            m = length(func());
        end
        
        function func_p = plusEpsShift(obj,i,vec_old,func)
            vec_p = vec_old;
            vec_p(i) = vec_p(i) + obj.eps;
            obj.vec = vec_p;
            func_p = func();
        end
        
        function func_n = minusEpsShift(obj,vec_old,func)
            vec_n = vec_old;
            vec_n(i) = vec_n(i) - obj.eps;
            obj.vec = vec_n;
            func_n = func();
        end
    end
    
end

