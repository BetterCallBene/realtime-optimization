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
        
        horizon;
        cDyn;
        ricM;
        
        s;
        lambda;
        q;
        mu;
    end
    
    methods
        
        function o = RealtimeSolver(varargin)
            % REALTIMESOLVER
            
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
            elseif(nargin == 6)
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
                
                o.lambda = varargin{3};
                o.s = varargin{4};
                o.q = varargin{5};
                o.mu = varargin{6};
                
            else
                error('wrong number of inputs');
            end
        end
        
        function [ res, est_y ] = fminrt(o, getLD, getLDD, n_timepoints)
            % FMINRT Solves the optimization problem in the SQP-Riccati approach
            
            %Initialize variables
            o.est_y = cell(n_timepoints,4);
            o.res = cell(n_timepoints,4);
            
            for i = 1:n_timepoints
                % Set estimated values
                vec = o.cDyn.getVecFromCells(o.s,o.q);
                o.cCost.vec = vec;
                o.cConst.vec = vec;
                
                %Update active set
                o.mu = o.cConst.checkIfActive(o.mu);
                
                % Calculate deltas with Riccati
                o.calculateSolution(i, getLD, getLDD);
                
                %Perform Newton Step and setup for next iteration
                o.performNewtonAndShift();
                
                %Estimate values at last timestep, by duplicating the values from the previous last step
                o.estimateNewHorizonPoint();
                
                %Save the estimated values
                o.storeEstimatedValues(i);
                
                %Display the actual timepoint in the command window
                disp(int2str(i));
            end
            
            res = o.res;
            est_y = o.est_y;
        end
    end
    
    methods(Test)
        
%         function testActiveSetChecker(o)
%             
%             hori = 20;
%             n_timepoints = 15;
%             o.setupTest(hori);
%             
%             %Initialize variables
%             o.est_y = cell(n_timepoints,4);
%             o.res = cell(n_timepoints,4);
%             
%             for i = 1:n_timepoints
%                 
%                 %Set values
%                 vec = o.cDyn.getVecFromCells(o.s,o.q);
%                 o.cCost.vec = vec;
%                 o.cConst.vec = vec;
%                 
%                 %Update active set
%                 o.mu = o.cConst.checkIfActive(o.mu);
%                 
%                 %Check if active set is correct
%                 for j = 1: hori
%                     aS_j = o.cConst.getActiveSet(j);
%                     %calculate inequalities
%                     is_leq_0 = o.cConst.get_ineq_con_at_t(j) >= 0;
%                     %TODO: refine
%                     o.assertEqual(sum(aS_j - is_leq_0) , 0);
%                     % TODO: Fehlen hier noch diverse Checks???
%                     %                     subplot(3,1,1), plot(aS_j,'r+-');
%                     %                     subplot(3,1,2), plot(is_leq_0,'r+-');
%                     %                     subplot(3,1,3), plot(o.cConst.get_ineq_con_at_t(j),'r+-');
%                 end
%                 
%                 % Calculate deltas with Riccati
%                 cLagrange = Lagrange();
%                 get_LD = @(cRTSolver, t) cLagrange.getLD(cRTSolver,t);
%                 get_LDD = @(cRTSolver,t) cLagrange.getLDD_approx(cRTSolver, t) ;
% 
%                 o.calculateSolution(i,get_LD, get_LDD);
%                 
%                 %Perform Newton Step and setup for next iteration
%                 o.performNewtonAndShift();
%                 
%                 %Estimate values at last timestep, by duplicating the values from the previous last step
%                 o.estimateNewHorizonPoint();
%                 
%                 %Save the estimated values
%                 o.storeEstimatedValues(i);
%                 
%                 %Display the actual timepoint in the command window
%                 disp(int2str(i));
%                 
%             end
%             
%         end
%         
        function testCalculateSolution(o)
            % TESTCALCULATESOLUTION Checks if the method calculateSolution
            % works correctly, using the \ operator and a tolerance of
            % 1e-15.
            
            % Initialize classes
            hori = 15;
            o.setupTest(hori);
            
            %Set values
            vec = o.cDyn.getVecFromCells(o.s,o.q);
            o.cCost.vec = vec;
            o.cConst.vec = vec;
            
            %Update active set
            o.mu = o.cConst.checkIfActive(o.mu);
            
            %Choose how to calculate the LDD (approximation or not)?
            cLagrange = Lagrange();
            get_LD = @(cRTSolver, t) cLagrange.getLD(cRTSolver,t);
            get_LDD = @(cRTSolver,t) cLagrange.getLDD_approx(cRTSolver, t) ;

            i = 1;
            tic;
            o.calculateSolution(i,get_LD, get_LDD);
            toc;
            
            %Build up delta
            delta_ric = o.buildUpDelta_ric();
           
            %Build up gradient and hesse
            [grad_L, hesse_L] = o.buildUpGradHesse(get_LD, get_LDD);
            
            %Solve with \
            delta_matlab = hesse_L \ grad_L;
            
            %Assert that both solutions are the same
            o.assertLessThan( norm(delta_matlab - delta_ric) , 1e-4 );
        end
        
        function testPerformNewtonAndShift(o)
            
            hori = 15;
            n_timepoints = 1; %TODO später wieder hochsetzen
            
            o.setupTest(hori);
            
            for i = 1: n_timepoints
                
                %Set values
                vec = o.cDyn.getVecFromCells(o.s,o.q);
                o.cCost.vec = vec;
                o.cConst.vec = vec;
                
                %Update active set
                o.mu = o.cConst.checkIfActive(o.mu);
                
                %Save old values
                old_s = o.s;
                old_lambda = o.lambda;
                old_q = o.q;
                old_mu = o.mu;
                
                %Choose how to calculate the LDD (approximation or not)?
                cLagrange = Lagrange();
                get_LD = @(cRTSolver, t) cLagrange.getLD(cRTSolver,t);
                get_LDD = @(cRTSolver,t) cLagrange.getLDD_approx(cRTSolver, t) ;
               %Build up delta and large gradient and hesse
                [grad_L, hesse_L] = o.buildUpGradHesse(get_LD, get_LDD);
                

                o.calculateSolution(i, get_LD, get_LDD);
                
                delta_ric = o.buildUpDelta_ric();
                % Perform Newton Step and setup for next iteration
                o.performNewtonAndShift();
              
                %Save the estimated values
                o.storeEstimatedValues(i);
                
                %Solve with \
                delta_matlab = hesse_L \ grad_L;
                
                % Compare results
                last = 0;
                n_var = 30  + sum(o.cConst.getActiveSet(1));
                o.assertLessThan( norm( old_lambda{1} + delta_matlab(1:13)- o.res{i,2}), 1e-4);
                o.assertLessThan( norm( old_s{1} + delta_matlab(14:26)- o.res{i,1}), 1e-4);
                o.assertLessThan( norm(old_q{1} + delta_matlab(27:30) - o.res{i,3}), 1e-4);
                if(n_var > 30)
                    delta_mu = zeros(8,1);
                    delta_mu(o.cConst.getActiveSet(1)) = delta_matlab(31:n_var);
                    o.assertLessThan( norm(old_mu( 1 : o.cConst.n_addConstr) + delta_mu - o.res{i,4}), 1e-4);
                end
                last = last + n_var;
                
                for l = 2:hori-1
                    n_var = 30  + sum(o.cConst.getActiveSet(l));
                    o.assertLessThan( norm(old_lambda{l} + delta_matlab(last +1 : last +13) - o.est_y{i,2}{l-1}), 1e-4 );
                    o.assertLessThan( norm(old_s{l} + delta_matlab(last+14 : last+26) - o.est_y{i,1}{l-1}), 1e-4);
                    o.assertLessThan( norm(old_q{l} + delta_matlab(last+27 : last+30) - o.est_y{i,3}{l-1}), 1e-4);
                    if(n_var > 30)
                        delta_mu = zeros(8,1);
                        delta_mu(o.cConst.getActiveSet(l-1)) = delta_matlab(last + 31:last +n_var);
                        o.assertLessThan( norm(old_mu( (l-1) * o.cConst.n_addConstr +1 : l * o.cConst.n_addConstr) + delta_mu  -o.est_y{i,4}( (l-2) *o.cConst.n_addConstr +1 : (l-1) * o.cConst.n_addConstr)), 1e-4);
                    end
                    last = last + n_var;
                end
                
            end
        end
    end
    methods(Access=private)
        function storeFirstIteration(o,i)
            o.res{i,1} = o.s{1} + o.ricM.delta_s{1};
            o.res{i,2} = o.lambda{1} + o.ricM.delta_lambda{1};
            o.res{i,3} = o.q{1} + o.ricM.delta_q{1};
            o.res{i,4} = o.mu(1:o.cConst.n_addConstr) + o.ricM.assembleMu(o.cConst.getActiveSet(1),1);
        end
        
        function storeEstimatedValues(o, i)
            o.est_y{i,1} = o.s;
            o.est_y{i,2} = o.lambda;
            o.est_y{i,3} = o.q;
            o.est_y{i,4} = o.mu;
        end
        
        function calculateSolution(o,i, getLD, getLDD)
            %Perform Riccati Steps
            for j = o.horizon+1:-1:1
                [LD, n_active_i] = getLD(o,j);
                LDD = getLDD(o,j);
                o.ricM.doStep(j,LDD, LD,n_active_i );
            end
            
            %Solve first Step
            o.ricM.solveStep(1);
            
            %Calculate new controls to pass it to the engines imidiatley
            %actualControl = q{1} + o.ricM.delta_q{1};
            
            %Store the result of the first iteration in res
            o.storeFirstIteration(i)
            
            %Solve the remaining steps to obtain a new iterate
            for j = 2:o.horizon+1
                o.ricM.solveStep(j);
            end
        end
        
        function performNewtonAndShift(o)
            % PERFORMNEWTONANDSHIFT Performs y + delta_y and shifts
            % every point to the left, as the actual timepoint increases by
            % one.
            for k = 2 : o.horizon
                o.s{k-1} = o.s{k} + o.ricM.delta_s{k};
                o.lambda{k-1} = o.lambda{k} + o.ricM.delta_lambda{k};
                o.q{k-1} = o.q{k} + o.ricM.delta_q{k};
                llI = (k-2) * o.cConst.n_addConstr +1 ;
                rlI = (k-1) * o.cConst.n_addConstr ;
                rrI = k * o.cConst.n_addConstr;
                o.mu(  llI : rlI ) = o.mu(rlI+1: rrI) +  o.ricM.assembleMu(o.cConst.getActiveSet(k) ,k);
            end
            k = o.horizon + 1 ;
            o.s{k-1} = o.s{k} + o.ricM.delta_s{k};
            o.lambda{k-1} = o.lambda{k} + o.ricM.delta_lambda{k};
        end
        
        function estimateNewHorizonPoint(o)
            % ESTIMATENEWHORIZONPOINT This method estimateds the the new
            % horizon point by duplicating the old horizon point.
            k = o.horizon +1 ;
            o.s{k} = o.s{k-1};
            o.lambda{k} = o.lambda{k-1};
            o.q{k-1} = o.q{k-2};
            o.mu((k-1) * o.cConst.n_addConstr +1 : end )= o.mu((k-2) * o.cConst.n_addConstr+1 : (k-1) * o.cConst.n_addConstr);
        end
        
        function [grad_L, hesse_L] = buildUpGradHesse(o,get_LD, get_LDD)
            hesse_L = zeros(1);
            grad_L= zeros(5,1);
            last = 0;
            
            for k = 1 : o.horizon
                n_var = 30  + sum(o.cConst.getActiveSet(k));
                hesse_L(last + 1 : last + 13 ,last  + 13 +1 :last  + 26 ) =  -eye(13);
                hesse_L(last  + 13 +1 : last  + 26 ,last  +1 :last + 13  ) =  -eye(13);
                hesse_L(last  +13+1 : last  +13 + n_var, last  +13 +1 : last +13 + n_var) = get_LDD(o,k);
                
                grad_L(last+1 : last + n_var) = get_LD(o, k);
                last = last + n_var;
            end
            k = o.horizon +1 ;
            
            LDk = get_LD(o, k);
            LDDk = get_LDD(o,k);
            
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
        
        function setupTest(o,hori)
            
            env = Environment();
            
            env.wind = @(s_t, t) s_t + 0.5 * [ones(3,1); zeros(10,1)];
            env.horizon = hori;
            env.setUniformMesh1(hori+1 ,1);
            cQ = Quadrocopter();
            
            
            tol = 1e-2;
            opts = odeset('RelTol',tol,'AbsTol',0.1*tol);
            cIntegrator = ode15sM(opts);
            cBQD = BasisQDyn(cQ,env,cIntegrator);
            
            o.cCost = CostsXU(cBQD,0.1, 50);
            multShoot = MultiShooting(cBQD);
            o.cConst = Constraints(multShoot);
            o.cDyn = o.cConst.dyn;
            o.horizon = o.cDyn.environment.horizon;
            o.ricM = RiccatiManager(o.horizon, o.cDyn.robot);
            
            o.s = cell(o.horizon +1,1);
            o.q = cell(o.horizon,1);
            o.lambda = cell(o.horizon +1 ,1);
            o.mu = ones( o.cConst.n_addConstr * (o.horizon+1),1);
            
            % Setup initial estimations
            for i = 1: o.horizon
                o.s{i} = [zeros(6,1); 1; zeros(6,1)];
                o.q{i} = zeros(4,1);
                o.lambda{i} = ones(cQ.n_state,1);
            end
            o.s{o.horizon +1} = [zeros(6,1); 1; zeros(6,1)];
            o.lambda{o.horizon +1} = ones(cQ.n_state,1);
        end
    end
end
