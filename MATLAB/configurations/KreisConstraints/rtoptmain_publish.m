currentpath= cd('..');
cd('..');
addpath(strcat(pwd, '/rtopt/'));
cd(currentpath);

PUBLISHABLE = false;

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
test = 20;

%% fmincon options    
% * Algorithm: interior-point
% * Display: iter
% * GradObj: on
% * GradConstr: on
% * Hessian: user-supplied
% * HessFcn: 
% * TolCon: $1e-06$ 
% * TolFun: $1e-06$
% * TolX: $1e-10$


n_int = 50;

%% Gitter und Intervallaenge
% * Intervallaenge: 50
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
% * States: $13$
% * Controls: $4$
% * Traegheitsmatrix: 
% * Gesamtmasse: $1.022$
% * kT: $1.5e-07$
% * kQ: $3e-09$
% * d: $0.22$
% * motor_m: $0.075$
% * motor_r: $0.015$


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
    %    \item Zeit: $80.0083$
    %    \item Exitflag: $-2$
    %    \begin{itemize}
    %       \item  1: First-order optimality measure was less than options
    %       TolFun and maximum constraint violation was less than options
    %       TolCon
    %       \item  0: Number of iterations exceeded options MaxIter or 
    %       number of function evaluations exceeded options MaxFunEvals.
    %       \item -2: No feasible point was found.
    %    \end{itemize}
    %    \item Iterations: $31$
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


 


