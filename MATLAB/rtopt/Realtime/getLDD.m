function [ H ] = getLDD(cCost, cConst, t )
% GETLDD 

if( t ~= cCost.dyn.environment.horizon +1 )
AB = cConst.get_eq_conD_block_t(t);
else
AB = [];
end
CD = cConst.get_ineq_conD_block_t(t);

ZERO = sparse(size(CD,1),size(CD,1));

%Active Set

activeSet_t = cConst.getActiveSet(t);

CD = CD(activeSet_t,:);
ZERO = ZERO(activeSet_t,activeSet_t);

costDall = cCost.get_costDD{1} ; %brauchen nur den Teil an der Stelle t
costD = costDall( (t-1)*(cConst.dyn.robot.n_state+cConst.dyn.robot.n_contr) + 1: t*(cConst.dyn.robot.n_state+cConst.dyn.robot.n_contr), (t-1)*(cConst.dyn.robot.n_state+cConst.dyn.robot.n_contr) + 1: t*(cConst.dyn.robot.n_state+cConst.dyn.robot.n_contr));

% cost_approx_t = cCost.get_costDD_approx(t);
% costD = cost_approx_t;

    
if (t ~= cCost.dyn.environment.horizon +1 )
H = [costD , CD' , AB' ;...
    CD , ZERO, sparse(size(ZERO,1),cConst.dyn.robot.n_state);...
    AB , sparse(cConst.dyn.robot.n_state,size(ZERO,1)) , sparse(cConst.dyn.robot.n_state,cConst.dyn.robot.n_state) ];
else
    H = [costD, CD'; CD, ZERO ];    
end


end

