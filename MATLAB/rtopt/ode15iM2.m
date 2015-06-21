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
        
        function [F, M, N, J] = helperCreateMatrizen(obj, y)
            
            [n_state, n_contr] = obj.getParams();
            
            n_constraints = 1;
            
            F = y(1:n_state, 1);
            M_vec = y(n_state + n_constraints +1 : n_state + n_constraints+ obj.M0_size, 1); %y(n_state + n_constraints +1 : n_state + n_constraints+ obj.M0_size, 1)
            N_vec = y(n_state + obj.M0_size + n_constraints + n_state + 1: n_state + obj.M0_size + n_constraints + n_state + obj.N0_size, 1);
            
            M = reshape(M_vec, [n_state, n_state]);
            N = reshape(N_vec, [n_state, n_contr]);
            
            J = [M, N];
        end
        
        function y0 = helperCreateInitialConditions(obj, varargin)
            [n_state, n_contr, n_var] = obj.getParams();
            
            if (nargin == 2)
                obj.nextStep(varargin{1})
            end
                        
            %y0 = obj.helperCreateVektor(obj.dyn.state(:, obj.timepoint), obj.M0, obj.N0);
            state = obj.vec_sav((obj.timepoint - 1) * n_var +1:(obj.timepoint - 1) * n_var + n_state);
            y0 = obj.helperCreateVektor(state, obj.M0, obj.N0);
        end
        
        function y = helperCreateVektor(obj, F, M, N)
            [n_state, n_contr, n_var] = obj.getParams();
            
            %y = zeros(n_state + n_constraints, 1);
            y = zeros(n_state, 1);
            y(1:n_state, 1) = F;
            %y(n_state+n_constraints:n_state+n_constraints, 1) = 0;
            %y(n_state + n_constraints +1 : n_state + n_constraints+ obj.M0_size, 1) = reshape(M, [obj.M0_size, 1]);
            %y(n_state + n_constraints+ obj.M0_size + 1: n_state + n_constraints+ obj.M0_size + n_state, 1) = sparse(n_state, 1);
            %y(n_state + obj.M0_size + n_constraints + n_state + 1: n_state + obj.M0_size + n_constraints + n_state + obj.N0_size, 1) = reshape(N, [obj.N0_size, 1]);
            %y(n_state + obj.M0_size + n_constraints + n_state + obj.N0_size + 1:end, 1) = sparse(n_contr, 1);
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
            
            [n_state] = obj.getParams();
            
            x = y(1:n_state);
            xp = yp(1:n_state);
            
            %qConstraint = 2* x(4:7)' * xp(4:7);
            %qConstraint = x(4)^2 + x(5)^2+x(6)^2 +x(7)^2 -1;
            obj.vec = [x; obj.u0];
                                                    %F,                     Jx * M
            f = obj.dyn.dot(obj.timepoint);
            
%             res = [ x(1:3) - f(1:3);
%                     x(4:7) .*xp(4:7) - f(4:7);
%                     x(8:end) - f(8:end)
%                   ];

            res = obj.mass(t, y)*xp - f;
            
        end
        
        function res = mass(obj, t, y)
            res = diag([ones(3, 1); y(4:7); ones(6, 1)]);
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
        
        
        function testOde(testCase)
            
            n_intervals = 2;
            timepoint = 1;
            
            testCase.n_intervalsInt = 50;
            
            [y0, yp0, old_interval, old_timepoint] = testCase.setupTest(n_intervals, timepoint);
            
            %JP = testCase.jpatternDy();
            %JPP = testCase.jpatternDyp();
            
            %testCase.opts = odeset('RelTol',1e-3,'AbsTol',1e-2,'Jacobian',@testCase.Jac,'Jpattern',{JP,JPP});
            %opts_ = odeset('RelTol',1e-4,'AbsTol',1e-3, 'Jacobian',@testCase.JacF);
            opts_ = odeset('RelTol',1e-4,'AbsTol',1e-5, 'Jacobian',@testCase.JacF);
            testCase.opts = opts_;
            
            tspan = [(timepoint -1)*testCase.h, timepoint*testCase.h];
            
            [y01,yp01] = decic(@testCase.funcToIntegrate,tspan(1),y0,[],yp0,[],opts_);
            testCase.odeTest(timepoint, y01, yp01);
            
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

