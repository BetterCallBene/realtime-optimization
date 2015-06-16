function [ H ] = getLDD(s_q,lambda,mu, cCost, cConst, t )
% GETLDD 


AB = cConst.get_eq_conD_block_t(t);
CD = cConst.get_ineq_conD_block_t(t);

%Active Set

activeSet_t = cConst.activeSet(


costDall = cCost.get_costDD ; %brauchen nur den Teil an der Stelle t
costD = costDall( (t-1)*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr) + 1: t*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr));
%wollte nicht zu viel in Cost umschreiben, daher such ich hier den
%passenden Wert raus
    

H = [costD , AB' ; AB , spzeros(cConstr.CountConstraints + cConstr.dyn.robot.n_state)];



end

