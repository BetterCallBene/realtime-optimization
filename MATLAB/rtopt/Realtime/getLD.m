function [ J ] = getLD( cCost, cConst, i )
%GETLD Calculates the 


% Find out size of the matrix



% Initialize J

constr = cConst.get_eq_con_at_t(i);
cost = cCost.get_costD;
cost = cost( (i-1)*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr) + 1: i*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr));


% Fill it up

J = [ cost ; constr( 2*cConstr.dyn.robot.n_state +1 : end) ; constr( cConstr.dyn.robot.n_state +1 :cConstr.dyn.robot.n_state) ];

end

