
%% Aufstellen des Problems




%TODO: options f�r fmincon festlegen
options                 = optimoptions('fmincon');
options.Algorithm       = 'interior-point';
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
    0, 0, 0,    ...     r           3      Ortsvektor
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
env.setUniformMesh(uint8(n_int));

cQ = Quadrocopter();

% TODO: Irgendwas besser f�r Startl�sung als rand finden
v0 = rand(cQ.n_var*(n_int+1),1);

% Initialisierung der Dynamik
%TODO: Klasse
cBQD = BasisQDyn(cQ, env);
CBQD.state = rand(cQ.n_state, n_int+1);
CBQD.contr = rand(cQ.n_contr, n_int+1);

% Wahl des Integrators
cFE = ForwEuler(cBQD);

% Initialisierung der Nebenbedingungen
%TODO: Klasse
cC = Constraints(cFE);

% Initialisierung Kostenfunktion
cCost = Costs(cBQD);

%Variablen f�r fmincon verf�gbar machen
%global objectCost objectConstr;
%objectCost = cCost;
%objectConstr = cC;

%% L�sung des Problems
tic;
%TODO: Realtime Ansatz einbauen, bzw. fmincon ersetzen/erweitern
v = fmincon(@costAdapter,v0,[],[],[],[],[],[],@constrAdapter, options);
toc;

%% Vorlage von Sebastian
% %% this file is used to call various subroutine for evaluating and testing
% %  during implementation 
% %  if evaluated in chronological order, all functions should provide 
% %  usable output
% 
% %%
% 
% options                 = optimoptions('fmincon');
% options.Algorithm       = 'interior-point';
% options.Display         = 'iter';
% options.GradObj         = 'on';
% options.GradConstr      = 'on';
% options.Hessian         = 'user-supplied';
% options.HessFcn         = @hessianAdapter;
% 
% %%
% n_int = 50;
% cDP = classDynParam([0.005;0.001], [2;1.5], [.4;.3]);
% xbc = [0.2, 0.5; -.2, .3; 0, 0; 0, 0];
% 
% cOCPp = classOCPparam(uint8(4),uint8(2),xbc,cDP);
% cOCPp.setUniformMesh(uint8(n_int));
% 
% %%
% v0 = rand(6*(n_int+1),1);
% 
% %%
% cR = classRobot();
% cR.param = cDP;
% 
% cD = classDyn(cR);
% 
% cFE = classForwEuler(cD,cOCPp);
% %cFE.state = rand(4,n_int+1);
% %cFE.contr = rand(2,n_int+1);
% %val = cFE.hDD();
% 
% cC = classConstraints(cFE,cOCPp);
% cC.vec = rand((n_int+1)*6,1);
% 
% [v1,v2,v3,v4] = cC.constr();
% 
% numDiffObj(cC);
% 
% cCost = classCosts(cOCPp);
% % cCost.vec = rand((n_int+1)*6,1);
% % cCost.get_costD();
% % numDiffObj(cCost);
% 
% %%
% global objectCost objectConstr;
% 
% objectCost = cCost;
% objectConstr = cC;
% 
% %%
% tic;
% v = fmincon(@costAdapter,v0,[],[],[],[],[],[],@constrAdapter, options);
% toc;



