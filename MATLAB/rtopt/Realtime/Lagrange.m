classdef Lagrange < TestEnv
    %LAGRANGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        cSolverRT; % For testing only
    end
    
    methods
        
        function L = getL(o,solverRT)
            %TODO: debug
            
            n_state = solverRT.cCost.dyn.robot.n_state;
            
            %Costs
            L = solverRT.cCost.get_cost();
            
            %Relaxed equality conditions
            %L = L + solverRT.lambda{1}' * (solverRT.cDyn.environment.wind(solverRT.s{1},1) - solverRT.s{1});
            for i = 1:solverRT.horizon+1
                tmp = solverRT.cConst.get_eq_con_at_t(i);
                L = L + solverRT.lambda{i}' * tmp( (i-1) * n_state +1 : i * n_state );
            end
            
            %Relaxed equality conditions
            for i = 1:solverRT.horizon
                activeSet_i = solverRT.cConst.getActiveSet(i);
                ineq_constr_i = solverRT.cConst.get_ineq_con_at_t(i);
                mu_i = solverRT.mu((i-1)*solverRT.cConst.n_addConstr +1: i * solverRT.cConst.n_addConstr);
                L = L + mu_i(activeSet_i)' * ineq_constr_i(activeSet_i);
            end
        end
        
        function [ J, n_active_i ] = getLD(o, solverRT, i )
            % GETLD Returns the NEGATIVE LD for timepoint i, because this
            % is the right hand side.
            % Initialize
            n_state = solverRT.cConst.dyn.robot.n_state;
            n_contr = solverRT.cConst.dyn.robot.n_contr;
            
            %Calculate eq and ineq constraints
            ineq_constr = solverRT.cConst.get_ineq_con_at_t(i);
            eq_constr = solverRT.cConst.get_eq_con_at_t(i);
            
            %Active Set (to find the right ineq constraints)
            activeSet_i = solverRT.cConst.getActiveSet(i);
            n_active_i = sum(activeSet_i);
            ineq_constr = ineq_constr(activeSet_i);
            
            %Calculate the cost derivative at the right place
            costDall = solverRT.cCost.get_costD;
            costD = costDall( (i-1)*(n_state+n_contr) + 1: i*(n_state+n_contr));
            
            %Berechne die Ableitung von den eqConstr nach s und q
            if i < (solverRT.horizon +1)
                eqConD = solverRT.cConst.get_eq_conD_block_t(i)' * solverRT.lambda{i+1};
                ineqConD = solverRT.cConst.get_ineq_conD_block_t(i);
                mu_i = solverRT.mu((i-1)*solverRT.cConst.n_addConstr +1: i *solverRT.cConst.n_addConstr);
                %   ineqConD = ineqConD(:, activeSet_i)' * mu_i(activeSet_i);
                ineqConD = ineqConD(activeSet_i,:)' * mu_i(activeSet_i);
            else
                eqConD = zeros(17,1);
                ineqConD = zeros(17,1);
            end
            
            if( i  > 1)
                lambda = [-solverRT.lambda{i};zeros(n_contr,1)];
                
            else
                lambda = zeros(n_state + n_contr,1);
            end
            % Fill it up
            J = -[ eq_constr( (i-1)*n_state +1 : i*n_state); costD + lambda + eqConD + ineqConD ; ineq_constr ];
            
            
        end
        
        function [ H ] = getLDD(o,solverRT, t )
            % GETLDD
            % Initialize
            n_state = solverRT.cConst.dyn.robot.n_state;
            n_contr = solverRT.cConst.dyn.robot.n_contr;
            
            % In the last block, there are no eq_conD considered
            if( t ~= solverRT.horizon +1 )
                AB = solverRT.cConst.get_eq_conD_block_t(t);
            else
                AB = [];
            end
            %Compute ineq_conD
            CD = solverRT.cConst.get_ineq_conD_block_t(t);
            ZERO = sparse(size(CD,1),size(CD,1)); %needed to fill matrix up
            
            %Active Set
            activeSet_t = solverRT.cConst.getActiveSet(t);
            CD = CD(activeSet_t,:);
            ZERO = ZERO(activeSet_t,activeSet_t);
            
            costDDall = solverRT.cCost.get_costDD{1} ; %brauchen nur den Teil an der Stelle t
            costDD = costDDall( (t-1)*(n_state+n_contr) + 1: t*(n_state+n_contr), (t-1)*(n_state+n_contr) + 1: t*(n_state+n_contr));
            
            %eq_conDDAll = solverRT.cConst.get_eq_conDD(solverRT.lambda{t+1});
            %eq_conDDAll = eq_conDDAll{1};
            eq_conDDAll = solverRT.cConst.get_eq_conDD_at_t(solverRT.lambda{t+1}, t);
            eq_conDD = eq_conDDAll( (t-1)*(n_state+n_contr) + 1: t*(n_state+n_contr), (t-1)*(n_state+n_contr)+ 1: t*(n_state+n_contr));
            
            % for approximation
            % cost_approx_t = solverRT.cCost.get_costDD_approx(t);
            % costD = cost_approx_t;
            
            
            if (t ~= solverRT.cCost.dyn.environment.horizon +1 )
                H = [costDD + eq_conDD, CD' , AB' ;...
                    CD , ZERO, sparse(size(ZERO,1),n_state);...
                    AB , sparse(n_state,size(ZERO,1)) , sparse(n_state,n_state) ];
            else
                H = [costDD + eq_conDD, CD'; CD, ZERO ];
            end
        end
        
        function [ H ] = getLDD_approx(o,solverRT, t )
            % GETLDD
            % Initialize
            n_state = solverRT.cConst.dyn.robot.n_state;
            n_contr = solverRT.cConst.dyn.robot.n_contr;
            
            % In the last block, there are no eq_conD considered
            if( t ~= solverRT.horizon +1 )
                AB = solverRT.cConst.get_eq_conD_block_t(t);
            else
                AB = [];
            end
            % Compute ineq_conD
            CD = solverRT.cConst.get_ineq_conD_block_t(t);
            ZERO = sparse(size(CD,1),size(CD,1)); %needed to fill matrix up
            
            % Active Set
            activeSet_t = solverRT.cConst.getActiveSet(t);
            CD = CD(activeSet_t,:);
            ZERO = ZERO(activeSet_t,activeSet_t);
            
            % Would be the exact solution
            %             costDall = solverRT.cCost.get_costDD{1} ; %brauchen nur den Teil an der Stelle t
            %             costD = costDall( (t-1)*(n_state+n_contr) + 1: t*(n_state+n_contr), (t-1)*(n_state+n_contr) + 1: t*(n_state+n_contr));
            
            cost_approx_t = solverRT.cCost.get_costDD_approx(t);
            costD = cost_approx_t;
            
            
            if (t ~= cCost.dyn.environment.horizon +1 )
                H = [costD , CD' , AB' ;...
                    CD , ZERO, sparse(size(ZERO,1),cConst.dyn.robot.n_state);...
                    AB , sparse(cConst.dyn.robot.n_state,size(ZERO,1)) , sparse(cConst.dyn.robot.n_state,cConst.dyn.robot.n_state) ];
            else
                H = [costD, CD'; CD, ZERO ];
            end
        end
        
    end
    
    
    methods(Test)
        
                function testgetLD(o)
                    horizon = 15;
                    o.setupTest(horizon);
        
                    func = @() o.getL(o.cSolverRT);
                    numDiff = o.numDiff_nD_AllT(func)';
        
                    for i = 1:horizon
        
                        anaDiff = o.getLD(o.cSolverRT,i);
                        %Perform checks
                        o.assertLessThan(norm(anaDiff + numDiff( (i-1) * 30 + 1 : i *30)), 1e-9);
        
                    end
                end
        
        function testgetLDD(o)
            horizon = 15;
            o.setupTest(horizon);
            
            for i =1:horizon
                
                func = @() o.getLD(o.cSolverRT,i)';
                anaDiff = o.getLDD(o.cSolverRT,i);
                numDiff1 = o.numDiff_nxnD(i,func);
                
                %Extend values, to derive after mu as well
                if i >= horizon
                   o.setupTest(horizon+1);
                   anaDiff = o.getLDD(o.cSolverRT,i);
                   numDiff1 = o.numDiff_nxnD(i,func);
                end
                
                numDiff2 = o.numDiff_nxnD(i+1,func);
                func2 = @() o.getLD(o.cSolverRT,i+1)';
                numDiff3 = o.numDiff_nxnD(i,func2);
                
                %Reshape the result
                numDiff2 = reshape(numDiff2(1,:,:), 30,30);
                numDiff3 = reshape(numDiff3(1,:,:),30,30);
                numDiff1 = reshape(numDiff1(1,:,:), 30,30);
                                
                %Perform checks
                QMR = anaDiff(1:17,1:17);
                AB = anaDiff(18:end, 1:17);
                ABtr = anaDiff(1:17, 18:end);
                
                o.assertLessThan( max(abs(anaDiff - anaDiff')), o.tol);
                o.assertLessThan( max(abs(QMR  + numDiff1(14:end, 14:end))),o.tol);
                if( i < horizon)
                    o.assertLessThan( max(abs(AB + numDiff2(1:13, 14:end))), o.tol);
                    o.assertLessThan( max(abs(ABtr + numDiff3(14:end, 1:13))), o.tol);
                end
                disp( [' Timestep: ', int2str(i), ' done'] );
            end
        end
        
        function testgetLD_withMus(o)
            %TODO: implement
        end
        
        function testgetLDD_withMus(o)
            %TODO: implement
        end
        
    end
    
    methods
        
        function setupTest(o,horizon)
            
            cEnv = Environment();
            cEnv.horizon = horizon;
            cEnv.wind = @(s_t ,t ) s_t + [1;1;1; zeros(10,1)];
            cEnv.setUniformMesh(uint16(horizon));
            cQ = Quadrocopter();
            cInt = ForwEuler();
            cBQD = BasisQDyn(cQ,cEnv,cInt);
            cMultShoot = MultiShooting(cBQD);
            cConst = Constraints(cMultShoot);
            cCost = CostsXU(cBQD,2,17);
            
            %Choose starting values
            s = cell(horizon +1,1);
            q = cell(horizon,1);
            lambda = cell(horizon +1 ,1);
            mu = ones( cConst.n_addConstr * (horizon+1),1);
            
            % Setup a initial estimation
            for i = 1: horizon
                s{i} = rand(13,1);
                q{i} = rand(4,1);
                lambda{i} = 10 * rand(cQ.n_state,1);
            end
            
            s{horizon +1} = rand(13,1);
            lambda{horizon +1} = 10 * rand(cQ.n_state,1);
            
            o.cSolverRT = RealtimeSolver(cCost, cConst,lambda, s, q, mu);
            
            vec = o.cSolverRT.cCost.dyn.getVecFromCells(s,q);
            o.cSolverRT.cCost.vec = vec;
            o.cSolverRT.cConst.vec = vec;
            
            o.cSolverRT.cConst.checkIfActive(mu);
            
        end
        
        function func_p = plusEpsShift(o, i,t,vec_old, func,n,dyn)
            vec_p = vec_old;
            if (i > 0 && i <=13)
                cell_index  = 1;
                j=i;
            elseif (i > 13 && i <= 26)
                cell_index = 2;
                j = i - 13;
                
            elseif (i > 26 && i <=30)
                cell_index = 3;
                j = i - 26;
                
            elseif  i > 30
                cell_index = 4;
            else
                error('i has infeasible value');
            end

            tmp = vec_p{cell_index};
            tmp2 = tmp{t};
            tmp2(j) = tmp2(j) + o.eps;
            tmp{t} = tmp2;
            vec_p{cell_index} = tmp;
            
            %Set value
            vec = dyn.getVecFromCells(vec_p{2}, vec_p{3});
            o.cSolverRT.cCost.vec = vec;
            o.cSolverRT.cConst.vec = vec;
            o.cSolverRT.lambda = vec_p{1};
            o.cSolverRT.s = vec_p{2};
            o.cSolverRT.q = vec_p{3};
            o.cSolverRT.mu = vec_p{4};
            
            % Evaluate function
            func_p = func();
            %Reset value
            vec = dyn.getVecFromCells(vec_old{2}, vec_old{3});
            o.cSolverRT.cCost.vec = vec;
            o.cSolverRT.cConst.vec = vec;
            o.cSolverRT.lambda = vec_old{1};
            o.cSolverRT.s = vec_old{2};
            o.cSolverRT.q = vec_old{3};
            o.cSolverRT.mu = vec_old{4};
        end
        
        function func_n =  minusEpsShift(o, i, t, vec_old, func, n, dyn)
            vec_n = vec_old;
            if (i > 0 && i <=13)
                cell_index  = 1;
                j=i;
            elseif (i > 13 && i <= 26)
                cell_index = 2;
                j = i - 13;
                
            elseif (i > 26 && i <=30)
                cell_index = 3;
                j = i - 26;
                
            elseif  i > 30
                cell_index = 4;
            else
                error('i has infeasible value');
            end
   
            tmp = vec_n{cell_index};
            tmp2 = tmp{t};
            tmp2(j) = tmp2(j) - o.eps;
            tmp{t} = tmp2;
            vec_n{cell_index} = tmp;
            
            %Set value
            vec = dyn.getVecFromCells(vec_n{2}, vec_n{3});
            o.cSolverRT.cCost.vec = vec;
            o.cSolverRT.cConst.vec = vec;
            o.cSolverRT.lambda = vec_n{1};
            o.cSolverRT.s = vec_n{2};
            o.cSolverRT.q = vec_n{3};
            o.cSolverRT.mu = vec_n{4};
            
            % Evaluate function
            func_n = func();
            %Reset value
            vec = dyn.getVecFromCells(vec_old{2}, vec_old{3});
            o.cSolverRT.cCost.vec = vec;
            o.cSolverRT.cConst.vec = vec;
            o.cSolverRT.lambda = vec_old{1};
            o.cSolverRT.s = vec_old{2};
            o.cSolverRT.q = vec_old{3};
            o.cSolverRT.mu = vec_old{4};
        end
        
        function [vec_old, n, m, n_timepoints, dyn] = setup(o, func, timepoint)
            vec_old = cell(4,1);
            vec_old{1} = o.cSolverRT.lambda;
            vec_old{2} = o.cSolverRT.s;
            vec_old{3} = o.cSolverRT.q;
            vec_old{4} = o.cSolverRT.mu;
            n_timepoints = o.cSolverRT.cDyn.environment.horizon;
            dyn = o.cSolverRT.cDyn;
            m = size(func(),1);
            n = 13 + 13+ 4;
        end
    end
end