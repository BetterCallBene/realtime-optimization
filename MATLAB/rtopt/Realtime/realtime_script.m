%% Define setting

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
horizon = 15;

env = Environment();
env.horizon = horizon;
env.xbc = xbc;

%Die Dynamik wird nur auf dem Horizon betrachtet:
env.setUniformMesh(uint16(horizon-1)); %TODO: Überprüfen, dass das nicht von den boundary conditions abhängt.

cQ = Quadrocopter();

% TODO: Irgendwas besser f�r Startl�sung als rand finden
%v0 = rand(cQ.n_var*horizon,1);

% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env);

% Wahl des Integrators
cFE = ForwEuler(cBQD);

% Initialisierung der Nebenbedingungen
cConst = Constraints(cFE);

% Initialisierung Kostenfunktion
cCost = Costs(cBQD);

%% Choose starting values

n_timepoints = 100 ; %How many timepoints, do we want to calculate.

s = cell(horizon +1,1);
q = cell(horizon,1);
lambda = cell(horizon +1 ,1);
mu = ones( cConst.n_addConstr * horizon,1);

% Setup a initial estimation
for i = 1: horizon 
s{i} = [zeros(6,1); 1; zeros(6,1)];
q{i} = zeros(4,1);
lambda{i} = ones(horizon,1);
end

s{horizon +1} = [zeros(6,1); 1; zeros(6,1)];
lambda{horizon +1} = ones(horizon,1);

%Choose how to calculate the LDD (approximation or not)?
getLDD = @(cost, const, t) getLDD(cost, const, t);

%% Calculate the solution with fminrt

tic;
%Use realtime solver
v = fminrt(cCost, cConst, getLDD, horizon, n_timepoints, s,q,lambda,mu);
toc;