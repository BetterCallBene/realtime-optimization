currentpath= cd('..');
cd('..');
addpath(strcat(pwd, '/rtopt/'));
cd(currentpath);

PUBLISHABLE = true;

global TEST;
TEST = false;

load('WindRand.mat', 'Wind');

%parpool(8);

options                 = optimoptions('fmincon');
options.Algorithm       = 'sqp';
options.Display         = 'iter';
options.GradObj         = 'on';
options.GradConstr      = 'on';
%options.Hessian         = 'user-supplied';
%options.HessFcn         = @hessianAdapter;
options.TolCon          = 1e-5;
options.TolFun          = 1e-5;
options.TolX            = 1e-6;
options.MaxFunEvals     = 60;
options.MaxIter         = 60;

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

xbc = [         ... Variablenname L�nge   Name
                ... Anfangsbedingung
                2, 0, 5,  ...     r           3      Ortsvektor
                1, 0, 0, 0, ...     q           4      Quaternion (Einheitsquaternion)
                0, 0, 0,    ...     v           3      Translatorische Geschwindigkeit
                0, 0, 0;    ...     w           3      Winkelgeschwindigkeit
                ... Endbedingung
                2, 0, 7,    ...
                1, 0, 0, 0, ...
                0, 0, 0,    ...
                0, 0, 0     ...
                ];


%Choose horizon
horizon = 49;
pointPerSecond = 1;

env = Environment();
env.horizon = horizon;
env.xbc = xbc;
env.wind = @(t, s_t , ctr) s_t;
%Die Dynamik wird nur auf dem Horizon betrachtet:
n_intervals = env.setUniformMesh1(horizon+1,pointPerSecond); 

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


% Wahl des Integrators
tol = 1e-2;
opts = odeset('RelTol',tol,'AbsTol',0.1*tol);
cIntegrator = ode15sM(opts); 
cIntegratorExt = ode15sM(opts);
% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env, cIntegrator);
cBQD.steadyPoint = [];

% Setup a initial estimation
steadyPoint = cBQD.steadyPoint;

% Initialisierung des Multiple Shootings
cMS = MultiShooting(cBQD);

% Initialisierung Kostenfunktion
%cCost = CostsComplet(cBQD, 1, 3, 0, 0);
cCost = CostsComplet(cBQD, 2.5, 75, 1, 1);

%Define Cam Position function
cCost.cam_pos = @(t) xbc(2, 1:3)'; %Standardwert verwenden [2, 0,
%5]

steadyPoint(1:3) = xbc(1, 1:3);
v0 = repmat(steadyPoint,(n_intervals+1), 1);
cBQD.vec =v0;  



% Initialisierung der Nebenbedingungen
cC = Constraints(cMS);

%Variablen fuer fmincon verfuegbar machen
global objectCost objectConstr;
objectCost = cCost;
objectConstr = cC;

%u_min = cQ.u_min;
%u_max = cQ.u_max;
%A = zeros(21, 17);
%A(14:17, 14:17) = -eye(4);
%A(18:end, 14:17) = eye(4);
%b = [zeros(13, 1);u_min *ones(4, 1); u_max *ones(4, 1)];

%repmat(A, 1, (n_intervals+1)); 

if PUBLISHABLE
    tic;
    [v, fval, exitflag, output] = fmincon(@costAdapter,v0,[],[],[],[],[],[],@constrAdapter, options);
    intervals = 1:(horizon +1);
    ProcessTime=toc;
    save('FData.mat', 'intervals', 'v', 'fval', 'exitflag', 'output', 'ProcessTime');
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
    load('FData.mat');
    Q = zeros(length(intervals), 12);
    array = reshape(v, [length(v)/length(intervals), length(intervals)])';

    Q(:, 1:3)  = array(:, 11:13); %omega
    [Q(:, 4), Q(:, 5), Q(:, 6)] = quat2angle((array(:, 4:7))); % Winkel
    Q(:, 7:9) =  array(:, 8:10); %v
    Q(:, 10:12) = array(:, 1:3); %position

    %% Beispiel: Norm of Quaternionen
    plot(intervals, array(:, 4).^2 + array(:, 5).^2+ array(:, 6).^2 + array(:, 7).^2);
end


 


