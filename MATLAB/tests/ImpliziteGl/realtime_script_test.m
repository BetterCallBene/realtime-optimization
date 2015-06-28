clear all
close all
%% Define setting

opts_ = odeset('RelTol',1e-5,'AbsTol',1e-7);

% Quadrocopter soll einen halben Meter nach unten fliegen
%TODO: Zusatzbedingungen, wie Anfangs und Endpunkt kommen in die Kostenfunktion
xbc = [         ... Variablenname L�nge   Name
                ... Anfangsbedingung
    0, 0, 0.5,  ...     r           3      Ortsvektor
    1, 0, 0, 0, ...     q           4      Quaternion (Einheitsquaternion)
    0, 0, 0,    ...     v           3      Translatorische Geschwindigkeit
    0, 0, 0;    ...     w           3      Winkelgeschwindigkeit
                ... Endbedingung
    0, 0, 0,    ...
    1, 0, 0, 0, ...
    0, 0, 0,    ...
    0, 0, 0     ...
];    

%Choose horizon
horizon = 10;

env = Environment();
env.horizon = horizon;
env.wind = @(s_t ,t ) s_t + [rand(3,1); zeros(10,1)];
env.xbc = xbc;

%Die Dynamik wird nur auf dem Horizon betrachtet:
env.setUniformMesh(uint16(horizon)); %TODO: Überprüfen, dass das nicht von den boundary conditions abhängt.

cQ = Quadrocopter();

% TODO: Irgendwas besser f�r Startl�sung als rand finden
v0 = rand(cQ.n_var*(horizon+1),1);

% Wahl des Integrators
integrator = ode15sM(opts_); %ode15iM2();%ForwEuler();


% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env,integrator);
cBQD.vec = v0;

% Initialisierung des Multiple Shootings
cMultShoot = MultiShooting(cBQD);

% Initialisierung der Nebenbedingungen
cConst = Constraints(cMultShoot);

% Initialisierung Kostenfunktion
cCost = Costs(cBQD);

%% Choose starting values

n_timepoints = 15 ; %How many timepoints, do we want to calculate.

s = cell(horizon +1,1);
q = cell(horizon,1);
lambda = cell(horizon +1 ,1);
mu = ones( cConst.n_addConstr * (horizon+1),1);

% Setup a initial estimation
for i = 1: horizon 
s{i} = [zeros(6,1); 1; zeros(6,1)];
q{i} = zeros(4,1);
lambda{i} = ones(cQ.n_state,1);
end

s{horizon +1} = [zeros(6,1); 1; zeros(6,1)];
lambda{horizon +1} = ones(cQ.n_state,1);

%Choose how to calculate the LDD (approximation or not)?
getLDD = @(cost, const, t) getLDD(cost, const, t);

% Initialisierung des Solvers
cRTSolver = RealtimeSolver(cCost, cConst);

%% Calculate the solution with fminrt

tic;
%Use realtime solver
[res, est_y  ] = cRTSolver.fminrt(getLDD, n_timepoints, s, q, lambda, mu);
toc;