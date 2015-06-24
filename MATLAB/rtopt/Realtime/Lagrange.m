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
            L = L + solverRT.lambda{1}' * solverRT.cDyn.environment.wind(solverRT.s{1},1);
            for i = 1:solverRT.horizon
                tmp = solverRT.cConst.get_eq_con_at_t(i);
                L = L + solverRT.lambda{i+1}' * tmp( (i-1) * n_state +1 : i * n_state );
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
            % Initialize J
            n_state = solverRT.cConst.dyn.robot.n_state;
            n_contr = solverRT.cConst.dyn.robot.n_contr;
            
            ineq_constr = solverRT.cConst.get_ineq_con_at_t(i);
            eq_constr = solverRT.cConst.get_eq_con_at_t(i);
            
            %Active Set
            
            activeSet_i = solverRT.cConst.getActiveSet(i);
            n_active_i = sum(activeSet_i);
            
            ineq_constr = ineq_constr(activeSet_i);
            
            
            costDall = solverRT.cCost.get_costD;
            costD = costDall( (i-1)*(n_state+n_contr) + 1: i*(n_state+n_contr));
            
            %berechne die Ableitung von den eqConstr nach s und q
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
            
            lambda = [-solverRT.lambda{i};zeros(n_contr,1)];
            
            % Fill it up
            J = -[ eq_constr( (i-1)*n_state +1 : i*n_state); costD + lambda + eqConD + ineqConD ; ineq_constr ];
            
        end
        
        function [ H ] = getLDD(o,solverRT, t )
            % GETLDD
            n_state = solverRT.cConst.dyn.robot.n_state;
            n_contr = solverRT.cConst.dyn.robot.n_contr;
            
            if( t ~= solverRT.horizon +1 )
                AB = solverRT.cConst.get_eq_conD_block_t(t);
            else
                AB = [];
            end
            CD = solverRT.cConst.get_ineq_conD_block_t(t);
            
            ZERO = sparse(size(CD,1),size(CD,1));
            
            %Active Set
            
            activeSet_t = solverRT.cConst.getActiveSet(t);
            
            CD = CD(activeSet_t,:);
            ZERO = ZERO(activeSet_t,activeSet_t);
            
            costDall = solverRT.cCost.get_costDD{1} ; %brauchen nur den Teil an der Stelle t
            costD = costDall( (t-1)*(n_state+n_contr) + 1: t*(n_state+n_contr), (t-1)*(n_state+n_contr) + 1: t*(n_state+n_contr));
            
            % cost_approx_t = solverRT.cCost.get_costDD_approx(t);
            % costD = cost_approx_t;
            
            
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
            %TODO: implement
            horizon = 12;
            o.setupTest(horizon);
            
            old_lambda = o.cSolverRT.lambda;
            old_s = o.cSolverRT.s;
            old_mu = o.cSolverRT.mu;
            old_q = o.cSolverRT.q;
            
            %             for i = 1:horizon+1
            i = 1;
            func = @() o.getL(o.cSolverRT);
            anaDiff = o.getLD(o.cSolverRT,i);
            numDiff = o.numDiff_nD_AllT(func);
            
            size(anaDiff)
            size(numDiff)
            
            
            %             end
            
        end
        
        function testgetLDD(o)
            %TODO implement
             horizon = 12;
            o.setupTest(horizon);
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
            cEnv.wind = @(s_t ,t ) s_t + [rand(3,1); zeros(10,1)];
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
                lambda{i} = rand(cQ.n_state,1);
            end
            
            s{horizon +1} = rand(13,1);
            lambda{horizon +1} = rand(cQ.n_state,1);
            
            o.cSolverRT = RealtimeSolver(cCost, cConst,lambda, s, q, mu);
            
            vec = o.cSolverRT.cCost.dyn.getVecFromCells(s,q);
            o.cSolverRT.cCost.vec = vec;
            o.cSolverRT.cConst.vec = vec;
            
            o.cSolverRT.cConst.checkIfActive(mu);
            
        end
        
        function func_p = plusEpsShift(o, i,t,vec_old, func,n,dyn)
            vec_p = vec_old;
            switch i
                case i > 0 && i <=13
                    cell_index  = 1;
                    j=i;
                case  i >13 && i <= 26
                    cell_index = 2;
                    j = i - 13;
                    
                case i > 26 && i <=30
                    cell_index = 3;
                    j = i - 26;
                    
                case  i > 30
                    cell_index = 4;
%                     
%                 otherwise
%                     error('i has infeasible value');
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
            switch i
                case i > 0 && i <=13
                    cell_index  = 1;
                    
                case  i >13 && i <= 26
                    cell_index = 2;
                    i = i - 13;
                    
                case i > 26 && i <=30
                    cell_index = 3;
                    i = i - 26;
                    
                case  i > 30
                    cell_index = 4;
                    
                otherwise
                    error('i has infeasible value');
            end
            tmp = vec_n{cell_index};
            tmp2 = tmp{t};
            tmp2(i) = tmp2(i) - o.eps;
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