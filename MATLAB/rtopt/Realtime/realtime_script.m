close all
%% Define setting

% Quadrocopter soll einen halben Meter nach unten fliegen

%Choose horizon
horizon = 12;
pointPerSecond = 1;


env = Environment();
env.horizon = horizon;
env.wind = @(s_t ,t ) s_t ; %+ [rand(3,1); zeros(10,1)];

%Die Dynamik wird nur auf dem Horizon betrachtet:
n_intervals = env.setUniformMesh1(horizon+1,pointPerSecond); 
cQ = Quadrocopter();

% Wahl des Integrators
tol = 1e-2;
opts = odeset('RelTol',tol,'AbsTol',0.1*tol);
%cIntegrator = ode15sM(opts);
cIntegrator = ForwEuler();

% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env,cIntegrator);

% Initialisierung des Multiple Shootings
cMultShoot = MultiShooting(cBQD);

% Initialisierung der Nebenbedingungen
cConst = Constraints(cMultShoot);

% Initialisierung Kostenfunktion
cCost = CostsXU(cBQD, 0.1, 50);

%% Choose starting values

n_timepoints = 15 ; %How many timepoints, do we want to calculate.

s = cell(n_intervals +1,1);
q = cell(n_intervals,1);
lambda = cell(n_intervals +1 ,1);
mu = ones( cConst.n_addConstr * (n_intervals+1),1);

% Setup a initial estimation
% TODO: Irgendwas besser f�r Startl�sung als rand finden

for i = 1: n_intervals 
s{i} = [zeros(3,1); 1; zeros(9,1)];
q{i} = 10000* ones(4,1);
lambda{i} = ones(cQ.n_state,1);
end

s{n_intervals +1} = [zeros(3,1); 1; zeros(9,1)];
lambda{n_intervals +1} = ones(cQ.n_state,1);

% Initialisierung des Solvers
cRTSolver = RealtimeSolver(cCost, cConst,lambda, s, q, mu);

%Choose how to calculate the LDD (approximation or not)?
cLagrange = Lagrange();
getLD = @(cRTSolver, t) cLagrange.getLD(cRTSolver,t);
getLDD = @(cRTSolver,t) cLagrange.getLDD(cRTSolver, t) ;

%% Calculate the solution with fminrt

tic;
%Use realtime solver
[res, est_y  ] = cRTSolver.fminrt(getLD, getLDD, n_timepoints);
toc;