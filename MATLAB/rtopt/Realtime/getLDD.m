function [ H ] = getLDD(s_q,lambda,mu, cCost, cConst, t )
% GETLDD 


AB = cConst.get_eq_conD_block_t(t);
CD = cConst.get_ineq_conD_block_t(t);

ZERO = sparse(size(CD,1),size(CD,1));

%Active Set

activeSet_t = cConst.getActiveSet(t);

CD = CD(activeSet_t,:);
ZERO = ZERO(activeSet_t,activeSet_t);

costDall = cCost.get_costDD ; %brauchen nur den Teil an der Stelle t
costD = costDall( (t-1)*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr) + 1: t*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr));
%wollte nicht zu viel in Cost umschreiben, daher such ich hier den
%passenden Wert raus
    

H = [costD , CD' , AB' ;...
    CD , ZERO, sparse(size(ZERO,1),cConstr.dyn.robot.n_state);...
    AB , sparse(cConstr.dyn.robot.n_state,size(ZERO,1)) , sparse(cConstr.dyn.robot.n_state,cConstr.dyn.robot.n_state) ];



end

