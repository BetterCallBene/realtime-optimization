clear
global TEST;
TEST = false;

%TODO: options f�r fmincon festlegen
options                 = optimoptions('fmincon');
options.Algorithm       = 'interior-point';
%options.Algorithm       = 'sqp';
options.Display         = 'iter';
options.GradObj         = 'on';
options.GradConstr      = 'on';
options.Hessian         = 'user-supplied';
options.HessFcn         = @hessianAdapter;

%TODO: n_int auf einen sinnvollen Wert festlegen
n_int = 50;
% Quadrocopter soll 5 Meter hoch fliegen
xbc = [         ... Variablenname L�nge   Name
                ... Anfangsbedingung
    0, 0, 0,  ...     r           3      Ortsvektor
    1, 0, 0, 0, ...     q           4      Quaternion (Einheitsquaternion)
    0, 0, 0,    ...     v           3      Translatorische Geschwindigkeit
    0, 0, 0;    ...     w           3      Winkelgeschwindigkeit
                ... Endbedingung
    0, 0, 5,    ...
    1, 0, 0, 0, ...
    0, 0, 0,    ...
    0, 0, 0     ...
];    

env = Environment();
env.xbc = xbc;
env.setUniformMesh(uint16(n_int));

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
cC = Constraints(cFE);

% Initialisierung Kostenfunktion
cCost = Costs(cBQD);

%Variablen f�r fmincon verf�gbar machen
global objectCost objectConstr;
objectCost = cCost;
objectConstr = cC;


%% L�sung des Problems
tic;
%TODO: Realtime Ansatz einbauen, bzw. fmincon ersetzen/erweitern
v = fmincon(@costAdapter,v0,A,b,[],[],[],[],@constrAdapter, options);
t = linspace(0, 1, n_int + 1);
save('../visualization/Circle.mat', 't', 'v')
toc;

