function [ineq_con,eq_con,ineq_conD,eq_conD] = constrAdapter(v)
% adapter for the fmincon interface for the constraints
% using the (global) object for constraints
global objectConstr;

% setting new input vector
objectConstr.vec = v;

tic;
[ineq_con, ineq_conD]     = objectConstr.get_ineq();

[eq_con, eq_conD]   = objectConstr.get_eq();

toc;