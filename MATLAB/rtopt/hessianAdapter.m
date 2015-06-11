function h = hessianAdapter(v,lambda)
% adapter for the fmincon interface for the hessian of the lagrangian
% using the (global) objects for cost and constraints

global objectCost objectConstr;

% setting new input vector
objectCost.vec = v;
objectConstr.vec = v;

% get Hessian for subfunctions
fDD     = objectCost.get_hess();
eqConDD   = objectConstr.get_eqhess();
inEqConDD = objectConstr.get_ineqhess();

% initialization and computation with langrage multipliers lambda
n = size(eqConDD,1);
m = size(inEqConDD, 1);

h = fDD{1};

for i=1:n
    h = h + lambda.eqnonlin(i)*eqConDD{i} ;
end

for i=1:m
    h = h + lambda.ineqnonlin(i)*inEqConDD{i};
end
