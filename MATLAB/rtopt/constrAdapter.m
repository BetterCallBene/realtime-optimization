function [ineq_con,eq_con,ineq_conD,eq_conD] = constrAdapter(v)
% adapter for the fmincon interface for the constraints
% using the (global) object for constraints
global objectConstr;

% setting new input vector
objectConstr.vec = v;

ineq_con    = [];
ineq_conD   = [];

eq_con   = objectConstr.get_func();
eq_conD  = objectConstr.get_jac();