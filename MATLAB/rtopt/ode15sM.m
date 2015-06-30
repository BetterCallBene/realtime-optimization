classdef ode15sM < Solver
    %FORWEULER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        y0;
    end
    methods
        
        function M15s = ode15sM(varargin)
            M15s@Solver(nargin, varargin); %Bug in Matlab
            M15s.opts.Jacobian = @M15s.Jac; 
        end
        
        function res = getJDot(obj)
            n_state = obj.dyn.robot.n_state;
            vec = obj.vec;
            state = vec(1:n_state);
            u0 = vec(n_state + 1 : end);
            res = obj.dyn.getJTilde(state, u0);
        end
        
        
        function y = integrate(obj, func, meshGrid, y0, yp0)
            %opts_ = odeset('RelTol',1e-4,'AbsTol',1e-6, 'InitialSlope', yp0, 'M15s.opts.InitialSlope );
            opts_ = obj.opts;
            opts_.InitialSlope = yp0;
            
            [t, y] = ode15s(func, meshGrid, y0, opts_);
            y = y(end, :)';
        end
        
        function [F, J, M, N] = ode(obj, timepoint, varargin)
            obj.nextStep(timepoint);
            
            y0_ = obj.helperCreateInitialConditions();
            yp0 = obj.helperCreateInitialConditionsDot();
            [F, J, M, N] = ode@Solver(obj, timepoint, y0_, yp0);
        end
        
        function dy = funcToIntegrate(obj, t, varargin)
            
            y = varargin{1};
            u0 = obj.u0;
            
            [state, M0, N0] = obj.helperCreateMatrizen(y);
            
            obj.vec = [state; u0];
            FTilde = obj.dyn.FTilde(y, u0);
            
            kM = obj.kM(M0); 
            kN = obj.kN(N0);
            
            dy = obj.helperCreateVektor(FTilde, kM, kN);
        end
        
        function yp = helperCreateInitialConditionsDot(obj)
            
            [n_state, n_contr, n_var] = getParams(obj);
            
            old_timepoint = obj.timepoint;
            obj.timepoint = 1;
            vec = obj.vec_sav((old_timepoint - 1) * n_var +1:(old_timepoint - 1) * n_var + n_var);
            obj.dyn.vec = vec;
            
            xp = obj.dyn.FTilde(vec(1:n_state, 1), vec(n_state+1:end));
            Mp = obj.kM(obj.M0);
            Np = obj.kN(obj.N0);
            
            obj.timepoint = old_timepoint;
            
            yp = obj.helperCreateVektor(xp, Mp, Np);
        end
        
        function dfdy = Jac(obj, t, y)
            
            u0 = obj.u0;
            [x0, M0, N0] = obj.helperCreateMatrizen(y);
            
            obj.vec = [x0; u0];
            
            %JTilde = obj.dyn.getJTilde(x0, u0);
            HTilde = obj.dyn.getHTilde(x0, u0);
            JTildex = obj.JDot_x; %JTilde(1:n_state, 1:n_state);
            
            dfdy = [obj.JacX(M0, N0, JTildex, HTilde), ...
                    obj.JacM(JTildex), ...
                    obj.JacN(JTildex) ...
                    ];
            
        end
        
        function dfdyX = JacX(obj, M0, N0, JTildex, HTilde)
            [n_state, n_contr] = obj.getParams();
            
            Mx = zeros(n_state * n_state, n_state);
            Nx = zeros(n_state * n_contr, n_state);
            
            %JTildeu = JTilde(1:n_state, n_state+1:end);
            
            for i = 1:n_state
                M0spalte = full(M0(:, i));
                HTildexx = HTilde(:, 1:n_state, 1:n_state);
                Mx((i-1)*n_state+1:i*n_state, :) = tprod(HTildexx, [1 -1 2], M0spalte, [-1]);
                if i <= n_contr
                    N0spalte = full(N0(:, i));
                    
                    Nx((i-1)*n_state+1:i*n_state, :) = tprod(HTildexx, [1 -1 2], N0spalte, [-1])  ...
                        + tprod(1, [-1],  HTilde(:, i + n_state, 1:n_state), [1 -1, 2]);
                end
            end
            
            dfdyX = [sparse(JTildex);
                sparse(Mx);
                sparse(Nx);
            ];
            
        end
        
        function dfdyM = JacM(obj, JTildex)
            Mdiag = JTildex;
            
            dfdyM = [
                    sparse(13, 169);
                    sparse(blkdiag(Mdiag, Mdiag, Mdiag, Mdiag, Mdiag,...
                Mdiag, Mdiag, Mdiag, Mdiag, Mdiag,....
                Mdiag, Mdiag, Mdiag ...
            ));
                sparse(4*13, 169);
            ];
            
        end
        
        function dfdyM = JacN(obj, JTildex)
            Ndiag = JTildex;
            
            dfdyM = [
                    sparse(13, 52);
                    sparse(169, 52);
                    sparse(blkdiag(Ndiag, Ndiag, Ndiag, Ndiag));
            ];
            
        end
        
        
        
        function [vec_old, n, m, timepoints, dyn]  = setup(obj,func)
            [n_state, n_contr, n_var] = getParams(obj);
            
            [old_interval] = obj.preToDo();
            
            %y0 = obj.helperCreateInitialConditions();
            
            y0_ = obj.y0;
            
            old_timepoint = obj.timepoint;
            %obj.u0 = obj.get_contr(old_timepoint);
            
            obj.nextStep(old_timepoint);
            
            timepoints = 0;
            dyn = obj.dyn;
            
            %n = 17;
            m = size(func([], y0_), 1);
            n = m;
            
            vec_old = y0_;
            
        end
        
        function  setupTest(obj,n_intervals, timepoint)
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
            %val = [zeros(3, 1); 1; 0 ;0 ;0; zeros(6, 1)];
            obj.u0 = [10000; 10000; 10000; 10000];
            %vec = [val;u];
            %obj.dyn.vec = [vec; vec; vec];
            
            %obj.dyn.vec = rand(17 * (n_intervals + 1), 1);
            
            obj.y0 = rand(234, 1);
            
        end
        
        %Overwrite function in TestEnv, because we don't want normed
        %quaternios here.
        function func_p = plusEpsShift(obj,i,t,vec_old,func, n_var, dyn)
            %vec_old = dyn.vec;
            vec_p = vec_old;
            t = 1;
            vec_p((t-1)* n_var + i) = vec_p((t-1)* n_var + i) + obj.eps;
            %dyn.backdoor_vec = vec_p;
            func_p = func([], vec_p);
            
        end
        
        %Overwrite function in TestEnv, because we don't want normed
        %quaternios here.
        function func_n = minusEpsShift(obj,i,t,vec_old,func, n_var, dyn)
            %vec_old = dyn.vec;
            vec_n = vec_old;
            t = 1;
            vec_n((t-1)* n_var + i) = vec_n((t-1)* n_var + i) - obj.eps;
            %dyn.backdoor_vec = vec_n;
            func_n = func([], vec_n);
        end
        
        
    end
    
    methods(Test)
        
        function testJac(testCase)
            n_intervals = 4;
            timepoint = 1;
            
            testCase.setupTest(n_intervals, timepoint);
            
            %testCase.timepoint = timepoint;
            numDiffJ = testCase.numDiff_nD(timepoint, @testCase.funcToIntegrate);
            %tic
            [anaJ] =testCase.Jac([], testCase.y0);
            %(anaJ(4:7, 4:7) - numDiffJ(4:7, 4:7))
            spy(abs(anaJ - numDiffJ)> 9e-3)
            %toc
            %testCase.assertLessThan(max(abs(anaJ - numDiffJ)),testCase.tol);
            %testCase.assertLessThan(max(abs(anaJD - numDiffJD)),testCase.tol);
        end
        
        
        
    end
    
    
    
    
end

