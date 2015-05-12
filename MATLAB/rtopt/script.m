
%% Aufstellen des Problems

%TODO: Klasse
env = classEnvConstant();
%TODO: Klasse
cR = classQuadrocopter(env);


%TODO: options für fmincon festlegen
options                 = optimoptions('fmincon');
options.Algorithm       = 'interior-point';
options.Display         = 'iter';
options.GradObj         = 'on';
options.GradConstr      = 'on';
options.Hessian         = 'user-supplied';
options.HessFcn         = @hessianAdapter;

%TODO: n_int auf einen sinnvollen Wert festlegen
n_int = 50;

%TODO: xbc=? Anfangs und Endbedingungen festlegen

%TODO: Klasse, cDP fliegt raus 
%cOCPp = classOCPparam(uint8(4),uint8(2),xbc,cDP);
cOCPp.setUniformMesh(uint8(n_int));

% TODO: Irgendwas besser für Startlösung als rand finden
v0 = rand(cR.n_state_contr*(n_int+1),1);

% Initialisierung der Dynamik
%TODO: Klasse
cD = classQuadrocopterDyn(cR);

% Wahl des Integrators
%TODO: Klasse anpassen
cFE = classForwEuler(cD,cOCPp);

% Initialisierung der Nebenbedingungen
%TODO: Klasse
cC = classConstraints(cFE,cOCPp);
%TODO: irgendwas besseres als rand finden
cC.vec = rand((n_int+1)*cR.n_state_contr,1);

% Initialisierung Kostenfunktion
cCost = classCosts(cOCPp);

%Variablen für fmincon verfügbar machen
global objectCost objectConstr;
objectCost = cCost;
objectConstr = cC;

%% Lösung des Problems
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



