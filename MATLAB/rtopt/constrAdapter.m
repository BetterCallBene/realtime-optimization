function [ineq_con,eq_con,ineq_conD,eq_conD] = constrAdapter(v)
% adapter for the fmincon interface for the constraints
% using the (global) object for constraints
global objectConstr;

% setting new input vector
objectConstr.vec = v;


%[ineq_con, ineq_conD]     = objectConstr.get_ineq();
ineq_con = [];
ineq_conD = [];

objectConstr.dode.flag_h = false; % bei jedem Schritt auf neue Hessematrix berechnen
[eq_con]   = objectConstr.get_eq_con();
[eq_conD]  = objectConstr.get_eq_conD();
