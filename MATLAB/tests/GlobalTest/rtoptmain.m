currentpath= cd('..');
cd('..');
addpath(strcat(pwd, '/rtopt/'));
cd(currentpath);

PUBLISHABLE = true;

global TEST;
TEST = false;

%parpool(8);

options                 = optimoptions('fmincon');
options.Algorithm       = 'interior-point';
options.Display         = 'iter';
options.GradObj         = 'on';
options.GradConstr      = 'on';
%options.Hessian         = 'user-supplied';
%options.HessFcn         = @hessianAdapter;
options.TolCon          = 1e-5;
options.TolFun          = 1e-6;
options.TolX            = 1e-5;
options.MaxFunEvals     = 1000000;
options.MaxIter         = 1000000;

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


%Choose horizon
horizon = 20;
pointPerSecond = 1;

env = Environment();
env.horizon = horizon;
env.wind = @(s_t ,t ) s_t + 0.1 * [rand(3,1); zeros(10,1)];
%Die Dynamik wird nur auf dem Horizon betrachtet:
n_intervals = env.setUniformMesh1(horizon+1,pointPerSecond); 


%% Environment
% Anfangsbedingung:

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
opts_ = odeset('RelTol',1e-2,'AbsTol',1e-3);
cIntegrator =  ode15sM(opts_); %ForwEuler(); %ode15sM(opts_); %ForwEuler(); %ode15sM(opts_); %% %ForwEuler();%ode15sM(opts_); %ForwEuler(); %ode15sM(opts_); %ode15sM(opts_);
% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env, cIntegrator);
cBQD.steadyPoint = [];
point = cBQD.steadyPoint;
% Setup a initial estimation
steadyPoint = cBQD.steadyPoint;

% Initialisierung Kostenfunktion
%cCost = CostsComplet(cBQD, 0.1, 2, 1, 1);
cCost = CostsXU(cBQD, 0.1, 2);
%Define Cam Position function
%cCost.cam_pos = @(t) cCost.skierCamPos(t);

steadyPoint(1:3) = cCost.cam_pos;
v0 = repmat(steadyPoint,(n_intervals+1), 1);
cBQD.vec =v0;   %rand(cQ.n_var * (n_int+1), 1);

cMS = MultiShooting(cBQD);

% Initialisierung der Nebenbedingungen
cC = Constraints(cMS);

%Variablen fuer fmincon verfuegbar machen
global objectCost objectConstr;
objectCost = cCost;
objectConstr = cC;

if PUBLISHABLE
    tic;
    [v, fval, exitflag, output] = fmincon(@costAdapter,v0,[],[],[],[],[],[],@constrAdapter, options);
    intervals = linspace(0, 1, n_intervals + 1);
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


 


