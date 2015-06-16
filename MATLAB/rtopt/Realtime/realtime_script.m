n_int = 50;
% Quadrocopter soll 5 Meter hoch fliegen
%TODO: Zusatzbedingungen kommen in die Kostenfunktion
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
env.setUniformMesh(horizon);

cQ = Quadrocopter();

% TODO: Irgendwas besser f�r Startl�sung als rand finden
v0 = rand(cQ.n_var*(n_int+1),1);


% Initialisierung der Dynamik
%TODO: Klasse
cBQD = BasisQDyn(cQ, env);
cBQD.vec = rand(cQ.n_var * (n_int+1), 1);

% Wahl des Integrators
cFE = ForwEuler(cBQD);

% Initialisierung der Nebenbedingungen
%TODO: Klasse
cConst = Constraints(cFE);

% Initialisierung Kostenfunktion
cCost = Costs(cBQD);

%% Choose starting value
%TODO: Hier muss man sich noch überlegen, wie man die zusätzlichen
%Constraints initialisiert...
y_var = 2*cQ.n_state + cQ.n_contr + size(cConst.get_eq_con(),1);
%TODO: init active set

y = zeros( y_var * horizon, 1);
y(7:y_var: horizon*y_var) = 1; % Normalize Quaternions
%%
tic;

%Use realtime solver
v = fminrt(cCost, cConst, horizon, y);
toc;