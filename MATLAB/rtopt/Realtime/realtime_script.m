close all
%% Define setting

% Quadrocopter soll einen halben Meter nach unten fliegen

%Choose horizon
horizon = 12;
pointPerSecond = 1;

env = Environment();
env.horizon = horizon;
%Die Dynamik wird nur auf dem Horizon betrachtet:
n_intervals = env.setUniformMesh1(horizon+1,pointPerSecond); 

cQ = Quadrocopter();

% Wahl des Integrators
tol = 1e-2;
opts = odeset('RelTol',tol,'AbsTol',0.1*tol);
cIntegrator = ode15sM(opts);
cIntegratorExt = ode15sM(opts);
%cIntegrator = ForwEuler();


cQExt = QuadrocopterExt(cQ, env, cIntegratorExt);
cQExt.steadyPoint = [];  %steadyPoint initialisieren: SteadyPoint ist eine globale Variable!!
cQExt.hForceExt = @(v) 0.1 * rand(3, 1) + cQ.getF_w(v);
cQExt.hMomentExt = @() 0.1 * rand(3, 1);
%Neue Windfunktion
env.wind = @(s_t, ctr)  cQExt.wind(s_t, ctr);
%env.wind = @(s_t ,t ) s_t + 0.1 * [rand(3,1); zeros(10,1)];
% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env, cIntegrator);

% Initialisierung des Multiple Shootings
cMultShoot = MultiShooting(cBQD);

% Initialisierung der Nebenbedingungen
cConst = Constraints(cMultShoot);

% Initialisierung Kostenfunktion
cCost = CostsComplet(cBQD, 0.1, 2, 1, 1);

n_timepoints = 20 ; %How many timepoints, do we want to calculate.

%Define Cam Position function
cCost.cam_pos = @(t) cCost.skierCamPos(t);

%% Choose starting values

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

% Initialisierung des Solvers
cRTSolver = RealtimeSolver(cCost, cConst,lambda, s, q, mu);

%Choose how to calculate the LDD (approximation or not)?
cLagrange = Lagrange();
getLD = @(cRTSolver, t) cLagrange.getLD(cRTSolver,t);
getLDD = @(cRTSolver,t) cLagrange.getLDD_approx_costDDpAlphaI(cRTSolver, t, 0.01 * ones(17,1) ) ;

%% Calculate the solution with fminrt

tic;
%Use realtime solver
[res, est_y  ] = cRTSolver.fminrt(getLD, getLDD, n_timepoints);
ProcessTime = toc;
