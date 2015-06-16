function [ res ] = fminrt(cCost, cConst, horizon, n_timepoints,s,q,lambda, mu)
% FMINRT Solves the optimization problem in the SQP-Riccati approach

%Check types of arguments
if ~isa(cConst,'Constraints')
    error('second variable must be of class Constraint');
end

%Initialize variables
cDyn = cConst.dyn;
ricM = RiccatiManager_woConstr(horizon, cDyn.robot);
res = cell(n_timepoints,4);

tic;
for i = 1:n_timepoints
    %TODO: Datentyp anpassen y ->  s_q, lambda_mu
    cCost.vec = y;
    cConst.vec = y;
    
    %TODO: Prüfe aktives Set von t - t+N-1
    mu = cConst.checkIfActive(mu);
    
    %Initialize RiccatiManager
    calcLDD = @(t) getLDD(s,q,lambda,mu,cCost, cConst,t);
    calcLD = @(t) getLD(y,cCost, cConst,t);
    
    
    %Perform Riccati Steps
    for j = horizon+1:-1:1
        ricM.doStep(j,calcLDD(j), calcLD(j));
    end
    
    %Solve first Steps
    ricM.solveStep(1);
    
    %Give new controls to the engines
    actualControl = q{1} + ricM.delta_q{1};
    
    %Solve the remaining steps to obtain a new iterate
    for j = 2:horizon+1
        ricM.solveStep(j);
    end
    
    %Perform Newton Step and setup for next iteration
    
    for k =2 : o.horizon
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
    
    %Estimate last timestep, by duplicating the previous last step
    s{k} = s{k-1};
    lambda{k} = lambda{k-1};
    q{k-1} = q{k-2};
    mu((k-1) * cConst.n_addConstr +1 : end} = mu((k-2) * cConst.n_addConstr+1 : (k-1) * cConst.n_addConstr);
    
    %Save result
    res(i,1) = s;
    res(i,2) = lambda;
    res(i,3) = q;
    res(i,4) = mu;
    
end
toc;
