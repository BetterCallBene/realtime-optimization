function [ H ] = getLDD(s_q,lambda,mu, cCost, cConst, t )
% GETLDD 


% Find out size of the matrix

% Initialize H



constrD = cConst.get_eq_conD_block_t(t);

costDall = cCost.get_costDD ; %brauchen nur den Teil an der Stelle t
costD = costDall( (t-1)*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr) + 1: t*(cConstr.dyn.robot.n_state+cConstr.dyn.robot.n_contr));
%wollte nicht zu viel in Cost umschreiben, daher such ich hier den
%passenden Wert raus
    

H = [costD , constrD' ; constrD , spzeros(cConstr.CountConstraints + cConstr.dyn.robot.n_state)];



end

