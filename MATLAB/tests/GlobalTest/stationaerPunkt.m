function stationaerPunkt()
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
%options.Hessian         = 'user-supplied';
%options.HessFcn         = @hessianAdapter;
options.TolCon          = 1e-4;
options.TolFun          = 1e-4;
options.TolX            = 1e-8;
options.MaxFunEvals        = 6000;
test = 20;

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


n_int = 20;

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
    0, 0, 0.5,    ...
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


% Wahl des Integrators
opts_ = odeset('RelTol',1e-5,'AbsTol',1e-6);
cIntegrator = ode15sM(opts_); %ode15sM(opts_);
% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env, cIntegrator);
vec = [zeros(3, 1); 1;0;0;0; zeros(6, 1);4087;4087;4087;4087];
cBQD.vec =repmat(vec, (n_int +1), 1); %rand(cQ.n_var * (n_int+1), 1); 

funcF = @(y) f(cBQD, y);
funcJ = @(y) J(cBQD, y);

options = optimoptions('fsolve','Display','iter', 'Algorithm', 'levenberg-marquardt'); % Option to display output
x = fsolve(funcF,vec, options);

end

function ret = f(cBQD, y)
    ret = cBQD.FTilde(y(1:13), y(14:end));
end

function ret = J(cBQD, y)
    ret = cBQD.getJTilde(y(1:13), y(14:end));
end

