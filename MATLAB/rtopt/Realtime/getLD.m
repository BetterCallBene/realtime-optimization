function [ J ] = getLD( cCost, cConst, i )



% Initialize J

constr = cConst.get_ineq_con_at_t(i);

%Active Set

activeSet_i = cConst.getActiveSet(i);

constr = constr(activeSet_i);


costall = cCost.get_costD;
cost = costall( (i-1)*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr) + 1: i*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr));


% Fill it up
if i == 1
    J = [ constr( 1 :cConstr.dyn.robot.n_state); cost ; constr ];
else
    J = [ constr( cConstr.dyn.robot.n_state +1 :2*cConstr.dyn.robot.n_state); cost ; constr ];
end


end

