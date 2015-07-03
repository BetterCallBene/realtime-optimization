classdef TestEnv < matlab.unittest.TestCase
    %TESTENV Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        eps = 1e-2;
        %tol = 1e-10;
        %Toleranz auf 1e-3 Heraufgesetzt, da x^3 Terme nur mit O(h^2)
        %angenaehert werden koennen
        tol = 1e-3;
        tolRK = 1e-5;
    end
    
    methods
        
        function numDiff = numDiff_nD(obj, timepoint, func)
            % NUMDIFFDND This method calculates numerically the derivative
            % of func at timepoint t, when func only depends on obj.vec and has m dim output
            [vec_old, n, m, timepoints, dyn] = obj.setup(func);
            numDiff = zeros(m,n);
            for i=1:n
                func_p = obj.plusEpsShift(i,timepoint,vec_old,func, n, dyn);
                func_n = obj.minusEpsShift(i,timepoint,vec_old,func, n, dyn);
                
                %Central difference
                numDiff(:,i) = (func_p - func_n)/2/obj.eps;
            end
        end
        
        function numDiff = numDiff_nD_AllT(obj, func)
            % NUMDIFF_NDALLT This method calculates numerically the
            % derivative of func, but for all timepoints
            
            [vec_old, n, m, n_timepoints, dyn] = obj.setup(func);
            numDiff = zeros(m , n * n_timepoints);
            
            for timepoint = 1:n_timepoints
                for i = 1:n
                    func_p = obj.plusEpsShift(i,timepoint,vec_old, func, n, dyn);
                    func_n = obj.minusEpsShift(i,timepoint,vec_old, func, n, dyn);
                    
                    %Central difference
                    numDiff( : , ((timepoint-1) * n) + i) ...
                        = (func_p - func_n)/2/obj.eps;
                    
                    
                end
                disp(int2str(timepoint));
            end
        end
        
        function numDiff = numDiff_nxnD(obj, timepoint, func)
            % NUMDIFF_DNXND This method calculates numerically the
            % derivative of func, when func only depends on obj.vec and has
            % n times n dimensional output
            [vec_old, n, m, n_timepoints, dyn] = obj.setup(func);
            numDiff = zeros(m,n,n);
            for i=1:n
                func_p = obj.plusEpsShift(i,timepoint, vec_old,func, n, dyn);
                func_n = obj.minusEpsShift(i,timepoint, vec_old,func, n, dyn);
                
                %Central difference
                numDiff(:,i,:) = (func_p - func_n )/2/obj.eps;
            end
        end
        
        function numDiff = numDiff_nxnD_AllT(obj, func)
            % NUMDIFF_NXND_ALLT This method calculates numerically the
            % derivative of func, when func only depends on obj.vec and has
            % n times n dimensional output
            [vec_old, n, m, n_timepoints, dyn] = obj.setup(func);
            
            numDiff  = zeros( m , n*n_timepoints,n*n_timepoints);
            for timepoint = 1:n_timepoints
                
                for i = 1:n
                    func_p = obj.plusEpsShift(i, timepoint, vec_old, func, n, dyn);
                    func_n = obj.minusEpsShift(i, timepoint, vec_old, func, n, dyn);
                    
                    %Central difference
                    numDiff(:,(timepoint -1 ) * n + i, : ) = (func_p - func_n)/2/obj.eps;
                end
            end
        end
        
        
        %Some help functions (typically overwritten in subclasses)
        function [vec_old, n, m, n_timepoints, dyn] = setup(obj, func, timepoint)
            vec_old = obj.vec;
            n_timepoints = obj.environment.n_timepoints;
            dyn = obj;
            n = obj.robot.n_var;
            m = size(func());
            m = m(1);
        end
        
        %Overwrite function in TestEnv, because we don't want normed
        %quaternios here.
        function func_p = plusEpsShift(obj,i,t,vec_old,func, n_var, dyn)
            %vec_old = dyn.vec;
            vec_p = vec_old;
            vec_p((t-1)* n_var + i) = vec_p((t-1)* n_var + i) + obj.eps;
            dyn.backdoor_vec = vec_p;
            func_p = func();
            dyn.backdoor_vec = vec_old;
        end
        
        %Overwrite function in TestEnv, because we don't want normed
        %quaternios here.
        function func_n = minusEpsShift(obj,i,t,vec_old,func, n_var, dyn)
            %vec_old = dyn.vec;
            vec_n = vec_old;
            vec_n((t-1)* n_var + i) = vec_n((t-1)* n_var + i) - obj.eps;
            dyn.backdoor_vec = vec_n;
            func_n = func();
            dyn.backdoor_vec = vec_old;
        end
    end
end

