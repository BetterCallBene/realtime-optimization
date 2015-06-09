function [ H ] = getLDD(s_q,lambda,mu, cCost, cConst, t )
% GETLDD 


% Find out size of the matrix

% Initialize H



constrD = cConst.get_eq_conD_block_t(t);

costD = cCost.get_eq.... ;
    

H = [costD , constrD' ; constrD , spzeros(cConstr.CountConstraints + cConstr.dyn.robot.n_state)];



end

