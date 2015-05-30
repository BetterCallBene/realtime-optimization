function [ res ] = fminrt( cCost, cConst, horizon,y,calcLD, calcLDD )
% FMINRT Solves the optimization problem in the SQP-Riccati approach

%Check types of arguments
if  ~isa(cCost, 'Costs')
    error('first variable must be of class Cost');
end
if ~isa(cConst,'Constraints')
    error('second variable must be of class Constraint');
end

%Initialize variables
cDyn = cCost.dyn;
n_lambda = size(cConst.get_eq_con(),1);
%n_mu = size(cConst.get_ineq_con(),1);
n_timepoints = cDyn.environment.n_timepoints;
res = zeros((cDyn.robot.n_var+n_lambda)*horizon, n_timepoints);

%Check if starting point y has correct length
if (length(y) / (cDyn.robot.n_var + n_lambda)) ~= horizon
    error('the starting point y has wrong length');
end

tic;
for i = 1:n_timepoints
   %Calculate the relevant matrices    
   LD = calcLD(i);
   LDD = calcLDD(i);
   
   %Solve LD + LDD*deltay = 0 with Riccati
   deltay = riccati2(LD, LDD, cDyn.environment.measurement(i));
   
   %Perform the Newtopn step
   y = y + deltay;
   
   %Save result
   res(:,i) = y;
   
   %Setup next iteration (for the estimation of the new horizon point we
   %duplicate the old horizon point)
   y = [y(cDyn.robot.n_var+n_lambda + 1 :end), y(end - (cDyn.robot.n_var + n_lambda)+1 , end)];
    
end
toc;
