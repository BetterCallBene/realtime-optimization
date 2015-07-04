clear all
close all
%% Setup
% Finds out how much time is between the moment the measurement is taken and
% the moment, when the control is calculated for the Riccati method and the
% MATLAB / operator.

horizon = 50;
pointPerSecond = 1;


% Select model
cQ = Quadrocopter();

% Initialize Environment with a little wind
env = Environment();
env.horizon = horizon;
env.wind = @(s_t ,t ) s_t + 0.1 * [ones(3,1); zeros(10,1)];
n_intervals = env.setUniformMesh1(horizon+1,pointPerSecond);

% Init Integrator
tol = 1e-4;
opts = odeset('RelTol',tol,'AbsTol',0.1*tol);
cIntegrator = ode15sM(opts);

% Init Dynamics
cBQD = BasisQDyn(cQ, env,cIntegrator);

% Initialisierung des Multiple Shootings
cMultShoot = MultiShooting(cBQD);

% Initialisierung der Nebenbedingungen
cConst = Constraints(cMultShoot);

% Initialisierung Kostenfunktion
cCost = CostsComplet(cBQD, 0.1, 2, 1, 1);
cCost.cam_pos = @(t) [2 ; 0 ;5];

s = cell(n_intervals +1,1);
q = cell(n_intervals,1);
lambda = cell(n_intervals +1 ,1);
mu = ones( cConst.n_addConstr * (n_intervals+1),1);

% Setup a initial estimation
steadyPoint = cBQD.steadyPoint;

% Place Quadrocopter at the desired camera position
steadyPoint(1:3) = cCost.cam_pos(1);

for i = 1: n_intervals
    s{i} = steadyPoint(1:cQ.n_state);
    q{i} = steadyPoint(cQ.n_state + 1 : cQ.n_var);
    lambda{i} = i *ones(cQ.n_state,1);
end

s{n_intervals +1} = 2 * steadyPoint(1:cQ.n_state);
lambda{n_intervals +1} = ones(cQ.n_state,1);

%Choose how to calculate the LDD (approximation or not)?
cLagrange = Lagrange();
getLD = @(cRTSolver, t) cLagrange.getLD(cRTSolver,t);
getLDD = @(cRTSolver,t) cLagrange.getLDD_approx_costDDpAlphaI(cRTSolver, t, zeros(17,1) ) ;

n_rep = 5;
timeRicTotal = zeros(n_rep,1);
timeBslTotal = zeros(n_rep,1);

timeRic_fromXtoQ = zeros(n_rep,1);
timeBsl_fromxtoQ = zeros(n_rep,1);



% Prepare measurement

vec = cBQD.getVecFromCells(s,q);
cCost.vec = vec;
cConst.vec = vec;
cCost.timepoint = 42;

mu = cConst.checkIfActive(mu);

ricM = RiccatiManager(horizon, cQ);

solverRT = RealtimeSolver(cCost, cConst, lambda, s, q, mu);


%% Get values for Riccati

% Execute code a several times, such that the variables are stored
% efficiently and that every field is initialized in RiccatiManager.


for i = 1:n_rep
    tot = tic;
    %Perform solution for one timepoint
    
    for j = horizon+1:-1:2
        [LD, n_active] = getLD(solverRT, j);
        LDD = getLDD(solverRT, j ) ;
        ricM.doStep(j,LDD, LD, n_active);
    end
    
    j = 1;
    
    fromXtoQ = tic;
    [LD, n_active] = getLD(solverRT, j);
    LDD = getLDD(solverRT, j ) ;
    ricM.doStep(j,LDD, LD, n_active);
    ricM.solveStep(1);
    act_q = q{1} + ricM.delta_q{1};
    timeRic_fromXtoQ(i) = toc(fromXtoQ);
    
    for j = 2:horizon + 1
        ricM.solveStep(j);
    end
    
    timeRicTotal(i) = toc(tot);
end


%% Get values for /

for i = 1:n_rep
    tot = tic;
    
    % Build up Hessian and Gradient
    hesse_L = 0;
    grad_L = 0;
    last = 0;
    
    for k = 1: horizon
        n_var = 30;
        if(k == 1)
            % Hesse_L can also be build from down to up, such that the measurement is needed in the last step.
            % This is too much effort, for a single comparison, so we are
            % adapting the time measurement.
            fromXtoQ = tic;
        end
        hesse_L(last + 1 : last + 13 ,last  + 13 +1 :last  + 26 ) =  -eye(13);
        hesse_L(last  + 13 +1 : last  + 26 ,last  +1 :last + 13  ) =  -eye(13);
        hesse_L(last  +13+1 : last  +13 + n_var, last  +13 +1 : last +13 + n_var) = getLDD(solverRT,k);
        
        grad_L(last+1 : last + n_var) = getLD(solverRT, k);
        last = last + n_var;
        
        if(k == 1)
            timeBsl_fromxtoQ(i) = toc(fromXtoQ);
        end
    end
    
    %Build up last part
    k = horizon + 1 ;
    LDk = getLD(solverRT, k);
    LDDk = getLDD(solverRT, k);
    
    hesse_L( last +1:last +13  , last + 13+1: last + 26 ) = -eye(13);
    hesse_L( last + 13+1: last + 26 , last +1: last +13) = -eye(13);
    hesse_L( last +13+1 : last + 13 + 13, last +13 +1 : last + 13 +13) = LDDk(1:13,1:13);
    
    grad_L(last +1 : last + 26 ) = LDk(1:26);
    
    
    % Solve the system
    fromXtoQ = tic;
    delta = hesse_L \ grad_L' ;
    % Compute the new control
    act_q = q{1} + delta(27:30);
    timeBsl_fromxtoQ(i) = timeBsl_fromxtoQ(i) + toc(fromXtoQ);
    
    
    timeBslTotal(i) = toc(tot);
end

