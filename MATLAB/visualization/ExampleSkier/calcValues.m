close all
%% Define setting

% Quadrocopter soll einen halben Meter nach unten fliegen

%Choose horizon
horizon = 15;
pointPerSecond = 1;

env = Environment();
env.horizon = horizon;
%Die Dynamik wird nur auf dem Horizon betrachtet:
n_intervals = env.setUniformMesh1(horizon+1,pointPerSecond); 

cQ = Quadrocopter();

% Wahl des Integrators
tol = 1e-2;
opts = odeset('RelTol',tol,'AbsTol',0.1*tol);
cIntegrator = ode15sM(opts);
cIntegratorExt = ode15sM(opts);
%cIntegrator = ForwEuler();


cQExt = QuadrocopterExt(cQ, env, cIntegratorExt);
cQExt.steadyPoint = [];  %steadyPoint initialisieren: SteadyPoint ist eine globale Variable!!
cQExt.hForceExt = @(v) 0.1 * rand(3, 1) + cQ.getF_w(v);
cQExt.hMomentExt = @() 0.1 * rand(3, 1);
%Neue Windfunktion
env.wind = @(t, s_t, ctr)  cQExt.wind(s_t, ctr);
%env.wind = @(s_t ,t ) s_t + 0.1 * [rand(3,1); zeros(10,1)];
% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env, cIntegrator);

% Initialisierung des Multiple Shootings
cMultShoot = MultiShooting(cBQD);

% Initialisierung der Nebenbedingungen
cConst = Constraints(cMultShoot);

% Initialisierung Kostenfunktion
cCost = CostsComplet(cBQD, 2.5, 2, 1, 1);

n_timepoints = 50 ; %How many timepoints, do we want to calculate.

%Define Cam Position function
cCost.cam_pos = @(t) cCost.skierCamPos_Short(t);

%% Choose starting values

s = cell(n_intervals +1,1);
q = cell(n_intervals,1);
lambda = cell(n_intervals +1 ,1);
mu = ones( cConst.n_addConstr * (n_intervals+1),1);

% Setup a initial estimation
steadyPoint = cBQD.steadyPoint;

% Place Quadrocopter at the desired camera position 
steadyPoint(1:3) = cCost.cam_pos(1);

for i = 1: n_intervals 
s{i} = steadyPoint(1:cQ.n_state);
q{i} = steadyPoint(cQ.n_state + 1 : cQ.n_var); 
lambda{i} = i *ones(cQ.n_state,1);
end

s{n_intervals +1} = 2 * steadyPoint(1:cQ.n_state);
lambda{n_intervals +1} = ones(cQ.n_state,1);

% Initialisierung des Solvers
cRTSolver = RealtimeSolver(cCost, cConst,lambda, s, q, mu);

%Choose how to calculate the LDD (approximation or not)?
cLagrange = Lagrange();
getLD = @(cRTSolver, t) cLagrange.getLD(cRTSolver,t);
getLDD = @(cRTSolver,t) cLagrange.getLDD_approx_costDDpAlphaI(cRTSolver, t, zeros(17,1) ) ;

%% Calculate the solution with fminrt

tic;
%Use realtime solver
[res, est_y  ] = cRTSolver.fminrt(getLD, getLDD, n_timepoints);
toc;
%%
cam_pos = zeros(3,1);
for i = 1:n_timepoints
    cam_pos(:,i) = cCost.cam_pos(i);
end


costF = 0;
pos = zeros(3,n_timepoints);
for i = 1:n_timepoints
    costF(i) = res{i,5};
    tmp = res{i,1};
    pos(:,i) = tmp(1:3);
end

norm_t = zeros(1,1200);
for i=1:n_timepoints
norm_t(i) = norm(cam_pos(:,i) - pos(:,i));
end
save('visualization/ExampleSkier/skierCamPos_Short.mat');

load('skierCamPos_Short.mat');


%% Plot functions
% Plot cost function
figure
plot(costF)
title('Cost function value');

% Plot CamPos and Pos
figure
plot3(cam_pos(1,:), cam_pos(2,:), cam_pos(3,:),'r')
hold on
plot3(pos(1,:), pos(2,:), pos(3,:),'b')
axis([-30 30  0 500 0 500])
view(-116,16);
title(' CamPosition (red) and computed position (blue)');


%Plot

figure
plot(norm_t);
title('Distance between camPosition and computed position');


%% For Screencast of downhill setting

pause(1);

hLarge = figure('name', 'Quadrocopter following a skier');
% Set size of window to fullHD
set(hLarge , 'Position', [0 0, 1920 , 1080 ]);


pause(5);

for i = 1:n_timepoints
    subplot(2,3,[1,2,4,5]);
   
    plot3(cam_pos(1,1:i), cam_pos(2,1:i), cam_pos(3,1:i),'b')
    hold on
    plot3(pos(1,1:i), pos(2,1:i), pos(3,1:i),'r');
    plot3(cam_pos(1,i), cam_pos(2,i), cam_pos(3,i),'b.', 'markersize', 7)
    plot3(pos(1,i), pos(2,i), pos(3,i),'r.', 'markersize', 15);
    axis([-20 20  0 400 0 400]);
    view(-116,16);
    title( 'Position of the drone(red) and the camera position(blue) in R^3');
    hold off
    
    subplot(2,3,3);
    plot(norm_t(1:i),'r');
    hold on 
    plot( i, norm_t(i), 'r.', 'markersize', 30);
    axis( [0 n_timepoints 0 16] );
    title( 'Distance between actual position and cam position');
    hold off
    
    
    subplot(2,3,6);
    plot(costF(1:i),'r');
    hold on 
    plot(i,costF(i), 'r.', 'markersize', 30);
    axis( [0 n_timepoints 0 1500])
    title( 'Cost function');
    hold off
    
    
    pause(0.1);
end

%% Show animation in nice framework

pos = zeros(n_timepoints,3);
ang = zeros(n_timepoints,4);
dot_pos = zeros(n_timepoints,3);
dot_ang = zeros(n_timepoints,3);

for i = 1:n_timepoints
   tmp  = res{i,1};
    pos(i,:) = tmp(1:3)';
    ang(i,:) = tmp(4:7)';
    dot_pos(i,:) = tmp(8:10)';
    dot_ang(i,:) = tmp(11:13)';
end

angle = zeros(n_timepoints, 3);
[angle(:, 1), angle(:, 2), angle(:, 3)] = quat2angle((ang));
% Prepare variables
global visualization

visualization.yout = [ dot_ang , angle , dot_pos , pos];
visualization.tout = ones(n_timepoints,1);

QuadAnim4;