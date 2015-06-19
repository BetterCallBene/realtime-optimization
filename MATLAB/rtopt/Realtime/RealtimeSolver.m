classdef RealtimeSolver < TestEnv
    % REALTIMESOLVER Just a dummy class for the fminrt function, which
    % enables us to write a proper unittests for fminrt and to refactor the
    % function fminrt in smaller functions, without creating additional
    % files.
    
    properties
        cCost;
        cConst;
        
        est_y;      %A cell to store the estimated values at i at {i,:}
        res;        %A cekk to store the result at time i at {i,:}
    end
    
    properties
        horizon;
        cDyn;
        ricM;
    end
    
    methods
        
        function o = RealtimeSolver(varargin)
            % REALTIMESOLVER
            
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
            elseif(nargin == 2)
                %Check types of arguments
                if ~isa(varargin{1},'Costs')
                    error('second variable must be of class Costs');
                end
                if ~isa(varargin{2},'Constraints')
                    error('second variable must be of class Constraints');
                end
                
                o.cCost = varargin{1};
                o.cConst = varargin{2};
                o.cDyn = o.cConst.dyn;
                o.horizon = o.cDyn.environment.horizon;
                o.ricM = RiccatiManager(o.horizon, o.cDyn.robot);
                
            else
                error('wrong number of inputs');
            end
        end
        
        function [ res, est_y ] = fminrt(o, getLDD, n_timepoints,s,q,lambda, mu)
            % FMINRT Solves the optimization problem in the SQP-Riccati approach
            
            %Initialize variables
            o.est_y = cell(n_timepoints,4);
            o.res = cell(n_timepoints,4);
            
            for i = 1:n_timepoints
                % Set estimated values
                vec = o.cDyn.getVecFromCells(s,q);
                o.cCost.vec = vec;
                o.cConst.vec = vec;
                
                %Update active set
                mu = o.cConst.checkIfActive(mu);
                
                % Calculate deltas with Riccati
                o.calculateSolution(i, getLDD, s, lambda, q, mu);
                
                %Perform Newton Step and setup for next iteration
                [s, lambda, q, mu] = o.performNewtonAndShift(s, lambda, q, mu);
                
                %Estimate values at last timestep, by duplicating the values from the previous last step
                [s, lambda, q, mu] = o.estimateNewHorizonPoint(s, lambda, q, mu);
                
                %Save the estimated values
                o.storeEstimatedValues(i,s, lambda, q ,mu);
                
                %Display the actual timepoint in the command window
                disp(int2str(i));
            end
            
            res = o.res;
            est_y = o.est_y;
        end
    end
    
    methods(Test)
        
        function testActiveSetChecker(o)
            
            hori = 12;
            n_timepoints = 50;
            [s, lambda, q, mu] = o.setupTest(hori);
            
            %Initialize variables
            o.est_y = cell(n_timepoints,4);
            o.res = cell(n_timepoints,4);
            
            for i = 1:n_timepoints
                
                %Set values
                vec = o.cDyn.getVecFromCells(s,q);
                o.cCost.vec = vec;
                o.cConst.vec = vec;
                
                %Update active set
                mu = o.cConst.checkIfActive(mu);
                
                %Check if active set is correct
                for j = 1: hori
                    aS_j = o.cConst.getActiveSet(j);
                    %calculate inequalities
                    is_leq_0 = o.cConst.get_ineq_con_at_t(j) >= 0;
                    %TODO: refine
                    o.assertEqual(sum(aS_j - is_leq_0) , 0);
                    subplot(3,1,1), plot(aS_j,'r+-');
                    subplot(3,1,2), plot(is_leq_0,'r+-');
                    subplot(3,1,3), plot(o.cConst.get_ineq_con_at_t(j),'r+-');
                end
                
                % Calculate deltas with Riccati
                get_LDD = @(cost, const, t) getLDD(cost, const, t);
                o.calculateSolution(i, get_LDD, s, lambda, q, mu);
                
                %Perform Newton Step and setup for next iteration
                [s, lambda, q, mu] = o.performNewtonAndShift(s, lambda, q, mu);
                
                %Estimate values at last timestep, by duplicating the values from the previous last step
                [s, lambda, q, mu] = o.estimateNewHorizonPoint(s, lambda, q, mu);
                
                %Save the estimated values
                o.storeEstimatedValues(i,s, lambda, q ,mu);
                
                %Display the actual timepoint in the command window
                disp(int2str(i));
            end
            
        end
        
        function testCalculateSolution(o)
            % TESTCALCULATESOLUTION Checks if the method calculateSolution
            % works correctly, using the \ operator and a tolerance of
            % 1e-15.
            
            % Initialize classes
            hori = 30;
            [s, lambda, q, mu] = o.setupTest(hori);
            
            %Set values
            vec = o.cDyn.getVecFromCells(s,q);
            o.cCost.vec = vec;
            o.cConst.vec = vec;
            
            %Update active set
            mu = o.cConst.checkIfActive(mu);
            
            %Choose how to calculate the LDD (approximation or not)?
            get_LDD = @(cost, const, t) getLDD(cost, const, t);
            i = 1;
            o.calculateSolution(i,get_LDD,s, lambda,q,mu);
            
            %Build up delta and large gradient and hesse
            [grad_L, hesse_L] = o.buildUpGradHesse(get_LDD);
            delta_ric = o.buildUpDelta_ric();
            
            %Solve with \
            delta_matlab = hesse_L \ grad_L;
            
            %Assert that both solutions are the same
            o.assertLessThan( norm(delta_matlab - delta_ric) , 1e-15 );
        end
    end
    
    methods(Access=private)
        function storeFirstIteration(o,i, s, lambda, q, mu)
            o.res{i,1} = s{1} + o.ricM.delta_s{1};
            o.res{i,2} = lambda{1} + o.ricM.delta_lambda{1};
            o.res{i,3} = q{1} + o.ricM.delta_q{1};
            o.res{i,4} = mu(1:o.cConst.n_addConstr) + o.ricM.assembleMu(o.cConst.getActiveSet(1),1);
        end
        
        function storeEstimatedValues(o, i, s, lambda, q, mu)
            o.est_y{i,1} = s;
            o.est_y{i,2} = lambda;
            o.est_y{i,3} = q;
            o.est_y{i,4} = mu;
        end
        
        function calculateSolution(o,i, getLDD,s,lambda,q,mu)
            %Perform Riccati Steps
            for j = o.horizon+1:-1:1
                [LD, n_active_i] = getLD(o.cCost,o.cConst,j);
                LDD = getLDD(o.cCost,o.cConst,j);
                o.ricM.doStep(j,LDD, LD,n_active_i );
            end
            
            %Solve first Step
            o.ricM.solveStep(1);
            
            %Calculate new controls to pass it to the engines imidiatley
            %actualControl = q{1} + o.ricM.delta_q{1};
            
            %Store the result of the first iteration in res
            o.storeFirstIteration(i, s, lambda, q, mu)
            
            %Solve the remaining steps to obtain a new iterate
            for j = 2:o.horizon+1
                o.ricM.solveStep(j);
            end
        end
        
        function [s, lambda, q, mu] = performNewtonAndShift(o, s, lambda, q, mu)
            % PERFORMNEWTONANDSHIFT Performs y + delta_y and shifts
            % every point to the left, as the actual timepoint increases by
            % one.
            for k = 2 : o.horizon
                s{k-1} = s{k} + o.ricM.delta_s{k};
                lambda{k-1} = lambda{k} + o.ricM.delta_lambda{k};
                q{k-1} = q{k} + o.ricM.delta_q{k};
                llI = (k-2) * o.cConst.n_addConstr +1 ;
                rlI = (k-1) * o.cConst.n_addConstr ;
                rrI = k * o.cConst.n_addConstr;
                mu(  llI : rlI ) = mu(rlI+1: rrI) +  o.ricM.assembleMu(o.cConst.getActiveSet(k) ,k);
            end
            k = o.horizon + 1 ;
            s{k-1} = s{k} + o.ricM.delta_s{k};
            lambda{k-1} = lambda{k} + o.ricM.delta_lambda{k};
        end
        
        function [s, lambda, q, mu] = estimateNewHorizonPoint(o, s, lambda, q, mu)
            % ESTIMATENEWHORIZONPOINT This method estimateds the the new
            % horizon point by duplicating the old horizon point.
            k = o.horizon +1 ;
            s{k} = s{k-1};
            lambda{k} = lambda{k-1};
            q{k-1} = q{k-2};
            mu((k-1) * o.cConst.n_addConstr +1 : end )= mu((k-2) * o.cConst.n_addConstr+1 : (k-1) * o.cConst.n_addConstr);
        end
        
        function [grad_L, hesse_L] = buildUpGradHesse(o,get_LDD)
            hesse_L = zeros(1);
            grad_L= zeros(5,1);
            last = 0;
            
            for k = 1 : o.horizon
                n_var = 30  + sum(o.cConst.getActiveSet(k));
                hesse_L(last + 1 : last + 13 ,last  + 13 +1 :last  + 26 ) =  -eye(13);
                hesse_L(last  + 13 +1 : last  + 26 ,last  +1 :last + 13  ) =  -eye(13);
                hesse_L(last  +13+1 : last  +13 + n_var, last  +13 +1 : last +13 + n_var) = get_LDD(o.cCost, o.cConst,k);
                
                grad_L(last+1 : last + n_var) = getLD(o.cCost, o.cConst, k);
                last = last + n_var;
            end
            k = o.horizon +1 ;
            
            LDk = getLD(o.cCost, o.cConst, k);
            LDDk = get_LDD(o.cCost, o.cConst,k);
            
            hesse_L( last +1:last +13  , last + 13+1: last + 26 ) = -eye(13);
            hesse_L( last + 13+1: last + 26 , last +1: last +13) = -eye(13);
            hesse_L( last +13+1 : last + 13 + 13, last +13 +1 : last + 13 +13) = LDDk(1:13,1:13);
            
            grad_L(last +1 : last + 26 ) = LDk(1:26);
        end
        
        function delta_ric = buildUpDelta_ric(o)
            delta_ric = [];
            for k = 1 : o.horizon
                y =  [o.ricM.delta_lambda{k}; o.ricM.delta_s{k}; o.ricM.delta_q{k}; o.ricM.delta_mu{k}];
                delta_ric = [delta_ric ;y];
            end
            k = o.horizon +1 ;
            y = [o.ricM.delta_lambda{k} ; o.ricM.delta_s{k}];
            delta_ric = [delta_ric ; y];
        end
        
        function [s, lambda, q, mu] = setupTest(o,hori)
            
            env = Environment();
            
            env.wind = @(s_t, t) s_t + 0.5 * [ones(3,1); zeros(10,1)];
            env.horizon = hori;
            env.setUniformMesh(uint16(hori));
            cQ = Quadrocopter();
            cBQD = BasisQDyn(cQ,env);
            cFE = ForwEuler(cBQD);
            o.cCost = Costs(cBQD);
            o.cConst = Constraints(cFE);
            o.cDyn = o.cConst.dyn;
            o.horizon = o.cDyn.environment.horizon;
            o.ricM = RiccatiManager(o.horizon, o.cDyn.robot);
            
            s = cell(o.horizon +1,1);
            q = cell(o.horizon,1);
            lambda = cell(o.horizon +1 ,1);
            mu = ones( o.cConst.n_addConstr * (o.horizon+1),1);
            
            % Setup initial estimations
            for i = 1: o.horizon
                s{i} = [zeros(6,1); 1; zeros(6,1)];
                q{i} = zeros(4,1);
                lambda{i} = ones(cQ.n_state,1);
            end
            s{o.horizon +1} = [zeros(6,1); 1; zeros(6,1)];
            lambda{o.horizon +1} = ones(cQ.n_state,1);
        end
    end
end
