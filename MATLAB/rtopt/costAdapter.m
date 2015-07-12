function [f,fD] = costAdapter(v)
% adapter for the fmincon interface for the cost function
% using the (global) object for cost 
global objectCost;

% setting new input vector
objectCost.vec = v;


f   = objectCost.get_func();
fD  = objectCost.get_jac();