function [ res, est_y ] = fminrt(cCost, cConst, getLDD, horizon, n_timepoints,s,q,lambda, mu)
% FMINRT Solves the optimization problem in the SQP-Riccati approach

%Check types of arguments
if ~isa(cConst,'Constraints')
    error('second variable must be of class Constraint');
end

%Initialize variables
cDyn = cConst.dyn;
ricM = RiccatiManager_woConstr(horizon, cDyn.robot);
est_y = cell(n_timepoints,4);
res = cell(n_timepoints,1);

tic;
for i = 1:n_timepoints    
    % Set estimated values 
    vec = cDyn.getVecFromCells(s,q);
    cCost.vec = vec;
    cConst.vec = vec;
    
    %Update active set
    mu = cConst.checkIfActive(mu);
    
    %Define handles to calculate LD, LDD 
    %TODO: das kann man wahrscheinlich aus der For Schleife nach oben
    %ziehen. Einfach mal ausprobieren, wenn der Alg funktioniert.
    calcLDD = @(t) getLDD(cCost, cConst,t);
    calcLD  = @(t) getLD(cCost, cConst,t);    
    
    %Perform Riccati Steps
    for j = horizon+1:-1:1
        ricM.doStep(j,calcLDD(j), calcLD(j));
    end
    
    %Solve first Step
    ricM.solveStep(1);
    
    %Calculate new controls to pass it to the engines imidiatley
    actualControl = q{1} + ricM.delta_q{1};
    
    %Store the result of the first iteration in res
    res(i,1) = s{1} + ricM.delta_s{1};
    res(i,2) = lambda{1} + ricM.delta_lambda{1};
    res(i,3) = actualControl;
    res(i,4) = mu(1:cConst.n_addConstr) + ricM.assembleMu(cConst.getActiveSet(1),1);
    
    %Solve the remaining steps to obtain a new iterate
    for j = 2:horizon+1
        ricM.solveStep(j);
    end
    
    %Perform Newton Step and setup for next iteration    
    for k = 2 : o.horizon
        s{k-1} = s{k} + ricM.delta_s{k};
        lambda{k-1} = lambda{k} + ricM.delta_lambda{k};
        q{k-1} = q{k} + ricM.delta_q{k};
        mu{k-1} = mu{k} + ricM.assembleMu(cConst.getActiveSet(k) ,k);
        llI = (k-2) * cConst.n_addConstr +1 ;
        rlI = (k-1) * cConst.n_addConstr ;
        rrI = k * cConst.n_addConstr;
        mu(  llI : rlI ) = mu(rlI+1: rrI) +  ricM.assembleMu(cConst.getActiveSet(k) ,k);
    end
    k = o.horizon + 1 ;
    s{k-1} = s{k} + ricM.delta_s{k};
    lambda{k-1} = lambda{k} + ricM.delta_lambda{k};
    
    %Estimate values at last timestep, by duplicating the values from the previous last step
    s{k} = s{k-1};
    lambda{k} = lambda{k-1};
    q{k-1} = q{k-2};
    mu((k-1) * cConst.n_addConstr +1 : end )= mu((k-2) * cConst.n_addConstr+1 : (k-1) * cConst.n_addConstr);
    
    %Save the estimated values
    est_y(i,1) = s;
    est_y(i,2) = lambda;
    est_y(i,3) = q;
    est_y(i,4) = mu;
    
end
toc;
%TODO: Some nice printout about running time, ...