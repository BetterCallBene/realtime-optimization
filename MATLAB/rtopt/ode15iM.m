classdef ode15iM < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        opts;
    end
    
    methods
        
        function odMi = ode15iM()
            %JP = JPat(N) | MvPat(N);?
            odMi@Solver();
            %obj.opts = odeset('RelTol',1e-5,'AbsTol',1e-4,'Jacobian',@odMi.Jac,'Jpattern',{JP,[]});
        end
        
        
        function yp = helperCreateInitialConditionsDot(obj)
            
            [n_state, n_contr, n_var] = getParams(obj);
            
            M0 = reshape(obj.M0, [obj.M0_size, 1]);
            N0 = reshape(obj.N0, [obj.N0_size, 1]);
            
            old_timepoint = obj.timepoint;
            obj.timepoint = 1;
            obj.dyn.vec =  obj.vec_sav((old_timepoint - 1) * n_var +1:(old_timepoint - 1) * n_var + n_var);
            
            xp = obj.dyn.dot(obj.timepoint);
            Mp = obj.kM(M0);
            Np = obj.kN(N0);
            
            obj.timepoint = old_timepoint;
            
            yp = obj.helperCreateVektor(xp, Mp, Np);
        end
        
        function y = integrate(obj, func, meshGrid, y0)
            yp0 = obj.helperCreateInitialConditionsDot();
            %% Sparsity pattern of df/dy
            
            y = ode15i(func, meshGrid, y0, yp0);
        end
        
        function res = funcToIntegrate(obj, t, varargin)
            y = varargin{1};
            yp = varargin{2};
            
            [n_state, n_contr] = obj.getParams();
            
            x = y(1:n_state);
            xp = yp(1:n_state);
            
            M_vec = y(n_state+1:n_state + obj.M0_size);
            N_vec = y(n_state + obj.M0_size+1:end);
            
            Mp_vec = yp(n_state+1:n_state + obj.M0_size);
            Np_vec = yp(n_state + obj.M0_size+1:end);
            
            M_ = reshape(M_vec, [n_state, n_state]);
            N_ = reshape(N_vec, [n_state, n_contr]);
            
            qConstraint = 2* x(4:7)' * xp(4:7);
            
            qdot_vec = zeros(1, n_state);
            qdot_vec(4:7) = 2.*xp(4:7);
            
            
            obj.vec = [x; obj.u0];
                                                    %F,                     Jx * M
            f = obj.dyn.dot(obj.timepoint);
            fM = reshape(obj.kM(M_vec), [obj.M0_size, 1]);%% Des da
            fN = reshape(obj.kN(N_vec), [obj.N0_size, 1]);
            
            qM = qdot_vec * M_;
            qN = qdot_vec * N_;
            
            res = [ xp - f;
                    qConstraint;
                    Mp_vec - fM;
                    qM';
                    Np_vec - fN;
                    qN';
                  ];
            
        end
        
        
        function [dfdy, dfdyp] = Jac(obj, t, y, yp)
            [n_state, n_contr] = obj.getParams();
            
            x = y(1:n_state);
            xp = yp(1:n_state);
            
            q_vec = zeros(1, n_state);
            qdot_vec = zeros(1, n_state);
            
            q_vec(4:7) = 2.*x(4:7);
            qdot_vec(4:7) = 2.*xp(4:7);
            
            M_vec = y(n_state+1:n_state + obj.M0_size);
            N_vec = y(n_state + obj.M0_size +1:end);
            M_ = reshape(M_vec, [n_state, n_state]);
            N_ = reshape(N_vec, [n_state, n_contr]);
            
            obj.vec = [x; obj.u0];
            
            JF = obj.dyn.dotD(obj.timepoint);
            HF = obj.dyn.dotDD(obj.timepoint);
            
            hugeM =zeros(n_state * n_state, n_state);
            hugeN =zeros(n_state * n_contr, n_state);
            
            k = 1;
            for i = 1:n_state
                for h = 1:n_state
                    hugeM(k, :) = (HF{h}(1:n_state, 1:n_state) * M_(:, i))';   
                    k = k + 1;
                end
            end
            k = 1;
            for i = 1:n_contr
                for h = 1:n_state
                    hugeN(k, :) = (HF{h}(1:n_state, 1:n_state) * N_(:, i))' +  HF{h}(n_state+i:n_state+i, 1:n_state);   
                    k = k + 1;
                end
            end
            
            JFx = JF(1:n_state, 1:n_state);
            bJFx1 = blkdiag(JFx, JFx, JFx, JFx, ...
                    JFx, JFx, JFx, JFx, ...
                    JFx, JFx, JFx, JFx, JFx);
            bJFx2 = blkdiag(JFx, JFx, JFx, JFx);
            
            s1 = obj.getS(n_state, xp); 
            s2 = obj.getS(n_contr, xp);
            
            spMQ1 = 2.* sparse(s1(:, 1), s1(:, 2), s1(:, 3), n_state, n_state * n_state);
            spMQ2 = 2.* sparse(s2(:, 1), s2(:, 2), s2(:, 3), n_contr, n_state * n_contr);
            
            qpM = [zeros(13, 3), ones( 13, 4), zeros( 13, 6)].* M_';
            qpN = [zeros( 4, 3), ones(4, 4), zeros( 4, 6)].* N_';
            
            dfdyX = [-JFx;
                qdot_vec;
                -hugeM;
                zeros(13, 13);
                -hugeN;
                zeros(4, 13);
                ]; 
            dfdyM = [zeros(14, 169);
                -bJFx1;
                spMQ1;
                zeros(56, 169)
                ];
            dfdyN = [zeros(196, 52);
                -bJFx2;
                spMQ2;
                ];
            dfdy = [dfdyX, dfdyM, dfdyN];
            
            dfdypX = [eye(13);
                q_vec;
                zeros(169, 13);
                2.*qpM;
                zeros(52, 13);
                2.*qpN
                ]; 
             dfdypM = [
                   zeros(14, 169);
                   eye(169);
                   zeros(69, 169);
                 ];
             dfdypN = [zeros(196, 52);
                 eye(52);
                 zeros(4, 52);
                 ];
            dfdyp = [dfdypX, dfdypM, dfdypN];
        end
        
    end
    methods
        
        function s = getS(obj, m, xp)
            s = zeros(4 * m, 3); %1->Zeile, 2->Spalte, 3->Wert
            s0 = 4;
            k = 1;
            for i = 1:m
                for j = 1:4
                    s(k, 1) = i;
                    s(k, 2) = s0;
                    s(k, 3) = xp(3+j);
                    s0 = s0 + 1;
                    k= k+1;
                end %-> 8 
                s0 = s0 + 13 -4;
            end
            
        end
        function res = jpatternDy(obj)
            
        end
    end
    methods  
        function [n, m] = setup1(obj,func, y0, yp0)
            
            n = size(y0, 1); %obj.dyn.robot.n_var;
            m = size(func([], y0, yp0), 1);
            
        end
        
        
        function func_p = plusEpsShift1(obj,i, y0, yp0, func)
            %vec_old = dyn.vec;
            vec_p = y0;
            vec_p(i) = vec_p(i) + obj.eps;
            func_p = func([], vec_p, yp0);            
        end
        
        function func_n = minusEpsShift1(obj,i, y0, yp0, func)
            %vec_old = dyn.vec;
            vec_p = y0;
            vec_p(i) = vec_p(i) - obj.eps;
            func_n = func([], vec_p, yp0);
        end
        
        function func_p = plusEpsShift2(obj,i, y0, yp0, func)
            %vec_old = dyn.vec;
            vec_p = yp0;
            vec_p(i) = vec_p(i) + obj.eps;
            func_p = func([], y0, vec_p);            
        end
        
        function func_n = minusEpsShift2(obj,i, y0, yp0, func)
            %vec_old = dyn.vec;
            vec_p = yp0;
            vec_p(i) = vec_p(i) - obj.eps;
            func_n = func([], y0, vec_p); 
        end
        
        function numDiff = numDiff_nD1(obj, timepoint, func, y0, yp0)
            % NUMDIFFDND This method calculates numerically the derivative
            % of func at timepoint t, when func only depends on obj.vec and has m dim output
            [n, m] = obj.setup1(func, y0, yp0);
            numDiff = zeros(m,n);
            for i=1:n
                func_p = obj.plusEpsShift1(i, y0, yp0, func);
                func_n = obj.minusEpsShift1(i, y0, yp0, func);
                
                %Central difference
                numDiff(:,i) = (func_p - func_n)/2/obj.eps;
            end
            
        end
        
        function numDiff = numDiff_nD2(obj, timepoint, func, y0, yp0)
            % NUMDIFFDND This method calculates numerically the derivative
            % of func at timepoint t, when func only depends on obj.vec and has m dim output
            [n, m] = obj.setup1(func, y0, yp0);
            numDiff = zeros(m,n);
            for i=1:n
                func_p = obj.plusEpsShift2(i, y0, yp0, func);
                func_n = obj.minusEpsShift2(i, y0, yp0, func);
                
                %Central difference
                numDiff(:,i) = (func_p - func_n)/2/obj.eps;
            end
            
        end
        
        function [y0, yp0, old_interval, old_timepoint] = setupTest(obj,n_intervals, timepoint)
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
            
            obj.dyn = BasisQDyn(model, env, obj);
            obj.dyn.vec = rand(17* (n_intervals+1), 1);
            
            obj.nextStep(timepoint);
            [old_interval] = obj.preToDo();
            
            y0 = obj.helperCreateInitialConditions();
            yp0 = obj.helperCreateInitialConditionsDot();
            
            old_timepoint = obj.timepoint;
            obj.u0 = obj.get_contr(old_timepoint);
            obj.timepoint = 1;
            
        end
        
        
    end
    
    methods(Test)
        
        function testJac(testCase)
            n_intervals = 4;
            timepoint = 3;
            
            [y0, yp0, old_interval, old_timepoint] = testCase.setupTest(n_intervals, timepoint);
            numDiffJ = testCase.numDiff_nD1(timepoint, @testCase.funcToIntegrate, y0, yp0);
            numDiffJD = testCase.numDiff_nD2(timepoint, @testCase.funcToIntegrate, y0, yp0);
            [anaJ, anaJD] = testCase.Jac([], y0, yp0);
            
            %spy(numDiffJ(:, 1:13) - anaJ(:, 1:13) > 1e-4)
%             figure
            %spy(abs(anaJ - numDiffJ) > 1e-4)
%             figure
%             spy(anaJ)
%             figure
%             spy(numDiffJ)
%             testCase.timepoint = old_timepoint;
%             testCase.postToDo(old_interval)
           
        end
        
        function testFuncToIntegrate(testCase)
            n_intervals = 3;
            timepoint = 2;
            
            testCase.setupTest(n_intervals);
            testCase.nextStep(timepoint);
            
            [old_interval] = testCase.preToDo();
            
            y0 = testCase.helperCreateInitialConditions();
            yp0 = testCase.helperCreateInitialConditionsDot();
            old_timepoint = testCase.timepoint;
            testCase.timepoint = 1;
            testCase.u0 = testCase.get_contr(old_timepoint);
            testCase.funcToIntegrate([], y0, yp0);
            
            testCase.timepoint = old_timepoint;
            
            testCase.postToDo(old_interval);
            
        end
        
        
        
    end
    
    
    
    
end

