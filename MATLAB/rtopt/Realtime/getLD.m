function [ J ] = getLD( cCost, cConst, i )
%GETLD Calculates the


% Find out size of the matrix



% Initialize J

constr = cConst.get_eq_con_at_t(i);
costall = cCost.get_costD;
cost = costall( (i-1)*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr) + 1: i*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr));


% Fill it up
if i == 1
    J = [ constr( 1 :cConstr.dyn.robot.n_state); cost ; constr( 2*cConstr.dyn.robot.n_state +1 : end) ];
else
    J = [ constr( cConstr.dyn.robot.n_state +1 :2*cConstr.dyn.robot.n_state); cost ; constr( 2*cConstr.dyn.robot.n_state +1 : end) ];
end


end

