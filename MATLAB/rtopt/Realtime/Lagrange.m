classdef Lagrange < TestEnv
    %LAGRANGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function [ J, n_active_i ] = getLD( solverRT, i )
            % Initialize J
            n_state = solverRT.cConst.dyn.robot.n_state;
            n_contr = solverRT.cConst.dyn.robot.n_contr;
            
            ineq_constr = solverRT.cConst.get_ineq_con_at_t(i);
            eq_constr = solverRT.cConst.get_eq_con_at_t(i);
            
            %Active Set
            
            activeSet_i = solverRT.cConst.getActiveSet(i);
            n_active_i = sum(activeSet_i);
            
            ineq_constr = ineq_constr(activeSet_i);
            
            
            costDall = cCost.get_costD;
            costD = costDall( (i-1)*(n_state+n_contr) + 1: i*(n_state+n_contr));
            
            %berechne die Ableitung von den eqConstr nach s und q
            if i < (solverRT.horizon +1)
                eqConD = solverRT.cConst.get_eq_conD_block_t(i)' * solverRT.lambda{i+1};
                ineqConD = solverRT.cConst.get_ineq_conD_block_t(i);
                mu_i = o.mu((i-1)*solverRT.cConst.n_addConstr +1: i *solverRT.cConst.n_addConstr);
                ineqConD = ineqConD(:, activeSet_i)' * mu_i(activeSet_i);
            else
                eqConD = zeros(17,1);
                ineqConD = zeros(17,1);
            end
            
            lambda = [-solverRT.lambda{i};zeros(n_contr,1)];
            
            % Fill it up
            J = -[ eq_constr( (i-1)*n_state +1 : i*n_state); costD + lambda + eqConD + ineqConD ; ineq_constr ];
            
        end
        
        function [ H ] = getLDD(solverRT, t )
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
    
end

