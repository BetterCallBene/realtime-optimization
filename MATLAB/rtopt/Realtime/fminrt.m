function [ res ] = fminrt(cCost, cConst, horizon,y)
% FMINRT Solves the optimization problem in the SQP-Riccati approach

%Check types of arguments
if ~isa(cConst,'Constraints')
    error('second variable must be of class Constraint');
end

%Initialize variables
cDyn = cConst.dyn;
n_lambda = size(cConst.get_eq_con(),1);
%n_mu = size(cConst.get_ineq_con(),1);
n_timepoints = cDyn.environment.n_timepoints;
res = zeros((cDyn.robot.n_var+n_lambda)*horizon, n_timepoints);

%Check if starting point y has correct length
if (length(y) / (cDyn.robot.n_var + n_lambda)) ~= horizon
    error('the starting point y has wrong length');
end

%TODO: aus dem y extrahieren
s;q;lambda; mu;
tic;
for i = 1:n_timepoints
    %TODO: Datentyp anpassen y ->  s_q, lambda_mu
    cCost.vec = y;
    cConst.vec = y;   
    
    %Initialize RiccatiManager
    ricM = RiccatiManager_woConstr(horizon, cDyn.robot);
    calcLDD = @(t) getLDD(s,q,lambda,mu,cCost, cConst,t);
    calcLD = @(t) getLD(y,cCost, cConst,t);
    
    %Perform Riccati Steps
    for j = horizon+1:-1:1
        ricM.doStep(j,calcLDD(j), calcLD(j));
    end
    
    %Solve first Steps
    ricM.solveStep(1);
    
    %Give new controls to the engines
    % doSomething
    
    %Solve the remaining steps to obtain a new iterate
    for j = 2:horizon+1
        ricM.solveStep(j);
    end
    
    %Perform Newton Step
    y = y + ricM.delta;
        
    %Save result
    res(:,i) = y;
    
    %Setup next iteration (for the estimation of the new horizon point we
    %duplicate the old horizon point)
    y = [y(cDyn.robot.n_var+n_lambda + 1 :end), y(end - (cDyn.robot.n_var + n_lambda)+1 , end)];
    
end
toc;
