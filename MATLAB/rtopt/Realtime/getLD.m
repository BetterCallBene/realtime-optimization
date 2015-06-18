function [ J, n_active_i ] = getLD( cCost, cConst, i )



% Initialize J

ineq_constr = cConst.get_ineq_con_at_t(i);
eq_constr = cConst.get_eq_con_at_t(i);

%Active Set

activeSet_i = cConst.getActiveSet(i);
n_active_i = sum(activeSet_i);

ineq_constr = ineq_constr(activeSet_i);


costall = cCost.get_costD;
cost = costall( (i-1)*(cConst.dyn.robot.n_state+cConst.dyn.robot.n_contr) + 1: i*(cConst.dyn.robot.n_state+cConst.dyn.robot.n_contr));


% Fill it up

%TODO: if weg
if i == 1
    J = -[ eq_constr( 1 :cConst.dyn.robot.n_state); cost ; ineq_constr ];
else
    %J = [ eq_constr( cConst.dyn.robot.n_state +1 :2*cConst.dyn.robot.n_state); cost ; ineq_constr ];
    
    J = -[ eq_constr( (i-1)*cConst.dyn.robot.n_state +1 : i*cConst.dyn.robot.n_state); cost ; ineq_constr ];
    
end


end

