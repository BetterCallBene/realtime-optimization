classdef ode15iM2 < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
        opts;
        dfdyPattern;
        dfdypPattern;
        
        dfdyPatternflag;
        dfdypPatternflag;
    end
    
    methods
        
        function odMi = ode15iM2()
            %JP = JPat(N) | MvPat(N);?
            odMi@Solver();
            odMi.dfdyPatternflag = true;
            odMi.dfdypPatternflag = true;
            
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
        
        function y = integrate(obj, func, meshGrid, varargin)
            y0 = varargin{1};
            if nargin == 4
                yp0 = obj.helperCreateInitialConditionsDot();
            else
                yp0 = varargin{2};
            end
            %% Sparsity pattern of df/dy
            opts_ = obj.opts;
            sol = ode15i(func, meshGrid, y0, yp0, opts_);
            y = sol.y(:, end);
        end
        
        function res = funcToIntegrate(obj, t, varargin)
            y = varargin{1};
            yp = varargin{2};
            
            u0 = obj.u0;
            [x0, M0, N0] = obj.helperCreateMatrizen(y);
            [xp0, M0Dot, N0Dot] = obj.helperCreateMatrizen( yp);
            
            [F, f, B] = obj.odeF(t, x0, xp0, u0);
            [M, N] = obj.odeMN(t, x0, u0, M0, M0Dot, N0, N0Dot, B, f);
            
            res = obj.helperCreateVektor(F, M, N);
        end
        
        
        
        function [res, f, B] = odeF(obj, t, x, xp, u0)
            
            obj.vec = [x; u0];
            
            B = obj.mass(t, x);
            f = obj.dyn.dot(obj.timepoint);
            
            res = B*xp - f;
        end
        
        function [M, N] = odeMN(obj, t, x0, u0, M0, M0dot, N0, N0dot, B, f)
            [n_state] = obj.getParams();
            
            Bdot = obj.massdot(t, x0);
            
            JTilde = obj.dyn.getJTilde(x0, u0);
            JTildex = JTilde(1:n_state, 1:n_state);
            JTildeu = JTilde(1:n_state, n_state+1:end);
            
            resM = zeros(13, 13);
            resN = zeros(13, 4);
            
            for i= 1:n_state
                row_iM = 0;
                row_iN = 0;
                for j = 1:n_state
                    row_iM = row_iM +  Bdot{i, j} * M0 * f(j) + B(i, j) * M0dot(j, :);
                    row_iN = row_iN +  Bdot{i, j} * N0 * f(j) + B(i, j) * N0dot(j, :);
                end
                resM(i, :) = row_iM  - JTildex(i, :) * M0;
                resN(i, :) = row_iN  - JTildex(i, :) * N0 - JTildeu(i, :);
            end
            M = reshape(sparse(resM), [obj.M0_size, 1]);
            N = reshape(sparse(resN), [obj.N0_size, 1]);
        end
        
        function res = mass(obj, t, y)
            res = diag([ones(3, 1); y(4:7); ones(6, 1)]);
        end
        
        function res = massdot(obj, t, y)
            res = cell(13, 13);
            for i = 1:13
                for j = 1:13
                    res{i, j} = sparse(1, 13);
                end
            end

            res{4, 4} = sparse(1, 4, 1, 1, 13);
            res{5, 5} = sparse(1, 5, 1, 1, 13);
            res{6, 6} = sparse(1, 6, 1, 1, 13);
            res{7, 7} = sparse(1, 7, 1, 1, 13);

        end
        
        function [dfdy, dfdyp] = JacF(obj, t, y, yp)
            [n_state, n_contr] = obj.getParams();
            
            x = y(1:n_state);
            %xp = yp(1:n_state);
            
            %q_vec = zeros(1, n_state);
            %qdot_vec = zeros(1, n_state);
            
            %q_vec(4:7) = 2.*x(4:7);
            %qdot_vec(4:7) = 2.*xp(4:7);

            obj.vec = [x; obj.u0];
            
            JF = obj.dyn.dotD(obj.timepoint);
            JFx = JF(1:n_state, 1:n_state);

            
            dfdy = -JFx;
            dfdyp = obj.mass(t,y);
        end
        
    end
    methods
        
%         function s = getS(obj, m, xp)
%             s = zeros(4 * m, 3); %1->Zeile, 2->Spalte, 3->Wert
%             s0 = 4;
%             k = 1;
%             for i = 1:m
%                 for j = 1:4
%                     s(k, 1) = i;
%                     s(k, 2) = s0;
%                     s(k, 3) = xp(3+j);
%                     s0 = s0 + 1;
%                     k= k+1;
%                 end %-> 8 
%                 s0 = s0 + 13 -4;
%             end
%             
%         end
        
%         function s = getSpattern(obj, m)
%             s = zeros(4 * m, 3); %1->Zeile, 2->Spalte, 3->Wert
%             s0 = 4;
%             k = 1;
%             for i = 1:m
%                 for j = 1:4
%                     s(k, 1) = i;
%                     s(k, 2) = s0;
%                     s(k, 3) = 1;
%                     s0 = s0 + 1;
%                     k= k+1;
%                 end %-> 8 
%                 s0 = s0 + 13 -4;
%             end
%             
%         end
        
%         function dfdy = jpatternDy(obj)
%             
%             if obj.dfdyPatternflag
%             
%                 [n_state, n_contr, n_var] = getParams(obj);
%                 
%                 JF = obj.dyn.dotDpattern();
%                 HF = obj.dyn.dotDDpattern();
%                 JFx = JF(1:n_state, 1:n_state);
% 
%                 bJFx1 = blkdiag(JFx, JFx, JFx, JFx, ...
%                         JFx, JFx, JFx, JFx, ...
%                         JFx, JFx, JFx, JFx, JFx);
%                 bJFx2 = blkdiag(JFx, JFx, JFx, JFx);
% 
%                 qdot_vec = [sparse(1, 3), ones(1, 4), sparse(1, 6)];
% 
%                 s1 = obj.getSpattern(n_state); 
%                 s2 = obj.getSpattern(n_contr);
% 
%                 spMQ1 = sparse(s1(:, 1), s1(:, 2), s1(:, 3), n_state, n_state * n_state);
%                 spMQ2 = sparse(s2(:, 1), s2(:, 2), s2(:, 3), n_contr, n_state * n_contr);
% 
% 
%                 M_ = ones(n_state, n_state);
%                 N_ = ones(n_state, n_contr);
% 
%                 hugeM =zeros(n_state * n_state, n_state);
%                 hugeN =zeros(n_state * n_contr, n_state);
% 
% 
%                 k = 1;
%                 for i = 1:n_state
%                     for h = 1:n_state
%                         hugeM(k, :) = (HF{h}(1:n_state, 1:n_state) * M_(:, i))';
%                         if i <= n_contr
%                             hugeN(k, :) = (HF{h}(1:n_state, 1:n_state) * N_(:, i))' +  HF{h}(n_state+i:n_state+i, 1:n_state); 
%                         end
%                         k = k + 1;
%                     end
%                 end
%                 spHugeM = sparse(hugeM);
%                 spHugeN = sparse(hugeN);
% 
% 
%                 dfdyX = [JFx;
%                     qdot_vec;
%                     spHugeM;
%                     sparse(13, 13);
%                     spHugeN;
%                     sparse(4, 13);
%                     ]; 
%                 dfdyM = [sparse(14, 169);
%                     bJFx1;
%                     spMQ1;
%                     sparse(56, 169)
%                     ];
%                 dfdyN = [sparse(196, 52);
%                     bJFx2;
%                     spMQ2;
%                     ];
%                 obj.dfdyPattern = logical([dfdyX, dfdyM, dfdyN]);
%                 obj.dfdyPatternflag = false;
%             end
%             dfdy = obj.dfdyPattern;
%         end
%         
%         function dfdyp = jpatternDyp(obj)
%             
%             if obj.dfdypPatternflag  
%                 q_vec = [sparse(1, 3), ones(1, 4), sparse(1, 6)];
%                 dfdypX = [speye(13);
%                     q_vec;
%                     sparse(169, 13);
%                     [sparse(13, 3), ones( 13, 4), sparse( 13, 6)];
%                     sparse(52, 13);
%                     [sparse( 4, 3), ones(4, 4), sparse( 4, 6)]
%                     ]; 
%                  dfdypM = [
%                        sparse(14, 169);
%                        speye(169);
%                        sparse(69, 169);
%                      ];
%                  dfdypN = [sparse(196, 52);
%                      speye(52);
%                      sparse(4, 52);
%                      ];
%                 obj.dfdypPattern = [dfdypX, dfdypM, dfdypN];
%                 obj.dfdypPatternflag = false;
%             end
%             dfdyp = obj.dfdypPattern;
%         end
    end
    methods  
        function [n, m] = setup1(obj,func, y0, yp0)
            
            n = size(y0, 1); %obj.dyn.robot.n_var;
            m = size(func([], y0, yp0), 1);
            
        end
        
        function res = FTildeTest(obj, t, y0, u0)
            res = obj.dyn.FTilde(y0, u0);
        end
        
        function res = JTildeTest(obj, t, y0, u0)
            res = obj.dyn.getJTilde(y0, u0);
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
            
            [n_state] = testCase.getParams();
            numDiffJ = testCase.numDiff_nD1(timepoint, @testCase.funcToIntegrate, y0, yp0);
            %numDiffJD = testCase.numDiff_nD2(timepoint, @testCase.funcToIntegrate, y0, yp0);
            
            %[anaJ, anaJD] = testCase.JacF([], y0, yp0);
            
            
            %testCase.assertLessThan(max(abs(anaJ - numDiffJ)),testCase.tol);
            %testCase.assertLessThan(max(abs(anaJD - numDiffJD)),testCase.tol);
 
            %testCase.timepoint = old_timepoint;
            %testCase.postToDo(old_interval);
            
            %patternDy = testCase.jpatternDy();
            %patternDyp = testCase.jpatternDyp();
            
            %patternDy - logical(anaJ)
            
        end
        
        function testJTilde(testCase)
            n_intervals = 4;
            timepoint = 3;
            [y0, yp0, old_interval, old_timepoint] = testCase.setupTest(n_intervals, timepoint);
            [n_state] = testCase.getParams();
            
            y01 = y0(1:13);
            u0 = testCase.get_contr(timepoint);
            numDiffJ = testCase.numDiff_nD1(timepoint, @testCase.FTildeTest, y01, u0);
            anaJ = testCase.JTildeTest([], y0, u0);
            
            testCase.assertLessThan(max(abs(anaJ(1:n_state, 1:n_state) - numDiffJ)),testCase.tol);
            
        end
        
        
        function testOde(testCase)
            
            n_intervals = 4;
            timepoint = 2;
            
            testCase.n_intervalsInt = 50;
            
            [y0, yp0, old_interval, old_timepoint] = testCase.setupTest(n_intervals, timepoint);
            [n_state] = testCase.getParams();
            
            opts_ = odeset('RelTol',1e-3,'AbsTol',1e-4);
            testCase.opts = opts_;
            
            tspan = [(timepoint -1)*testCase.h, timepoint*testCase.h];
            tic;
            [y01,yp01] = decic(@testCase.funcToIntegrate,tspan(1),y0,[],yp0,[],opts_);
            
            absSchaetzCalc= norm(y01 - y0, 1);
            disp('Abstand des geschaetzen zudem berechneten Wert');
            disp(absSchaetzCalc);
            
            [F, J] = testCase.odeTest(timepoint, y01, yp01);
            toc
            Q = norm(F(4:7));
            % Differenz zur 1
            testCase.assertLessThan(Q - 1,testCase.tolRK);
            
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

