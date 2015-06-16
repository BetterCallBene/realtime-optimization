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
            obj.opts = odeset('RelTol',1e-5,'AbsTol',1e-4,'Jacobian',@odMi.Jac,'Jpattern',{JP,[]});
        end
        
        
        function yp = helperCreateInitialConditionsDot(obj)
            M0 = reshape(obj.M0, obj.M0_size);
            N0 = reshape(obj.N0, obj.N0_size);
            
            obj.dyn.vec =  obj.vec;
            
            xp = obj.dyn.dot(obj.timepoint);
            Mp = obj.kM(M0);
            Np = obj.kN(N0);
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
            
            qdot_vec = zeros(1, n_state);
            qdot_vec(4:7) = 2.*xp(4:7);
            
            obj.vec = [x; obj.u0];
                                                    %F,                     Jx * M
            f = obj.dyn.dot(obj.timepoint);
            fM = reshape(obj.kM(M_vec), [obj.M0_size, 1]);%% Des da
            fN = reshape(obj.kN(N_vec), [obj.N0_size, 1]);
            
            M_ = reshape(Mp_vec, [n_state, n_state]);
            N_ = reshape(Np_vec, [n_state, n_contr]);
            
            qM = qdot_vec * M_;
            qN = qdot_vec * N_;
            
            res = [ xp - f;
                    Mp_vec - fM;
                    qM;
                    Np_vec - fN;
                    qN;
                  ];
            
        end
        
        function [dfdy, dfdyp] = Jac(obj, t, y, yp)
            [n_state] = obj.getParams();
            
            x = y(1:n_state);
            xp = yp(1:n_state);
            
            qdot_vec = zeros(1, n_state);
            qdot_vec(4:7) = 2.*xp(4:7);
            
            M_vec = y(n_state+1:n_state + obj.M0_size);
            M_ = reshape(M_vec, [n_state, n_state]);
            
            JF = obj.dyn.dotD(obj.timepoint);
            HF = obj.dyn.dotDD(obj.timepoint);
            
            obj.vec = [x; obj.u0];
            
            hugeM =zeros(n_state * n_state, n_state);
            
            for h = 1:n_state
                hugeM((h-1)*n_state + 1:(h-1)*n_state + n_state, :) = HF{h} * M_;                
            end
            
            dfdy = [-JF;
                qdot_vec;
                -hugeM;
                zeros(13, 13);
                
                ]; 
            dfdyp = [];
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

