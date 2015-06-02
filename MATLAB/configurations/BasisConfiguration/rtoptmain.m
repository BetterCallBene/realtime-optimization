currentpath= cd('..');
cd('..');
addpath(strcat(pwd, '/rtopt/'));
cd(currentpath);

PUBLISHABLE = true;

global TEST;
TEST = false;

options                 = optimoptions('fmincon');
options.Algorithm       = 'interior-point';
options.Display         = 'iter';
options.GradObj         = 'on';
options.GradConstr      = 'on';
options.Hessian         = 'user-supplied';
options.HessFcn         = @hessianAdapter;
options.TolCon          = 1e-6;
options.TolFun          = 1e-6;
options.TolX            = 1e-10;

%% fmincon options    
% * Algorithm: $(options.Algorithm)$
% * Display: $(options.Display)$
% * GradObj: $(options.GradObj)$
% * GradConstr: $(options.GradConstr)$
% * Hessian: $(options.Hessian)$
% * HessFcn: $(options.HessFcn)$
% * TolCon: $$(options.TolCon)$$ 
% * TolFun: $$(options.TolFun)$$
% * TolX: $$(options.TolX)$$



n_int = 50;

%% Gitter und Intervallaenge
% * Intervallaenge: $(n_int)$
%

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

%% Environment
% Anfangsbedingung:


env = Environment();
env.xbc = xbc;
env.setUniformMesh(uint16(n_int));
  
cQ = Quadrocopter();

%% Quadrocopter Eigenschaften    
% * States: $$(cQ.n_state)$$
% * Controls: $$(cQ.n_contr)$$
% * Traegheitsmatrix: 
% * Gesamtmasse: $$(cQ.m)$$
% * kT: $$(cQ.kT)$$
% * kQ: $$(cQ.kQ)$$
% * d: $$(cQ.d)$$
% * motor_m: $$(cQ.motor_m)$$
% * motor_r: $$(cQ.motor_r)$$


v0 = rand(cQ.n_var*(n_int+1),1);

% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env);
cBQD.vec = rand(cQ.n_var * (n_int+1), 1);

% Wahl des Integrators
cFE = ForwEuler(cBQD);

% Initialisierung der Nebenbedingungen
cC = Constraints(cFE);

% Initialisierung Kostenfunktion
cCost = Costs(cBQD);

%Variablen fuer fmincon verfuegbar machen
global objectCost objectConstr;
objectCost = cCost;
objectConstr = cC;

if PUBLISHABLE
    tic
    [v, fval, exitflag, output] = fmincon(@costAdapter,v0,[],[],[],[],[],[],@constrAdapter, options);
    intervals = linspace(0, 1, n_int + 1);
    ProcessTime=toc;
    save('Data.mat', 'intervals', 'v', 'fval', 'exitflag', 'output', 'ProcessTime');
    %% Ergebnis
    %
    % <latex>
    % \begin{itemize}
    %    \item Zeit: $$(ProcessTime)$$
    %    \item Exitflag: $$(exitflag)$$
    %    \begin{itemize}
    %       \item  1: First-order optimality measure was less than options
    %       TolFun and maximum constraint violation was less than options
    %       TolCon
    %       \item  0: Number of iterations exceeded options MaxIter or 
    %       number of function evaluations exceeded options MaxFunEvals.
    %       \item -2: No feasible point was found.
    %    \end{itemize}
    %    \item Iterations: $$(output.iterations)$$
    % \end{itemize}
    % </latex>
else
    load('Data.mat');
    Q = zeros(length(intervals), 12);
    array = reshape(v, [length(v)/length(intervals), length(intervals)])';

    Q(:, 1:3)  = array(:, 11:13); %omega
    [Q(:, 4), Q(:, 5), Q(:, 6)] = quat2angle((array(:, 4:7))); % Winkel
    Q(:, 7:9) =  array(:, 8:10); %v
    Q(:, 10:12) = array(:, 1:3); %position

    %% Beispiel: Norm of Quaternionen
    plot(intervals, array(:, 4).^2 + array(:, 5).^2+ array(:, 6).^2 + array(:, 7).^2);
end


 


