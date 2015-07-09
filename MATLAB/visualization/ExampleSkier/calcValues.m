%% Define setting
close all

% Quadrocopter soll einem Skifahrer folgen

%Choose horizon
horizon = 18;
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
env.wind = @(t, s_t, ctr)  cQExt.wind(t,s_t, ctr);
%env.wind = @(t, s_t, ctr ) s_t + 0.1 * [rand(3,1); zeros(10,1)];
% Initialisierung der Dynamik
cBQD = BasisQDyn(cQ, env, cIntegrator);

% Initialisierung des Multiple Shootings
cMultShoot = MultiShooting(cBQD);

% Initialisierung der Nebenbedingungen
cConst = Constraints(cMultShoot);

% Initialisierung Kostenfunktion

cCost = CostsComplet(cBQD, 1, 0.25, 1, 1);

n_timepoints = 4*60 ; %How many timepoints, do we want to calculate.

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


costF = zeros(n_timepoints,1);
pos = zeros(3,n_timepoints);
contr = zeros(4,n_timepoints);
norm_lambda = zeros(1,n_timepoints);
for i = 1:n_timepoints
    costF(i) = res{i,5};
    tmp = res{i,1};
    pos(:,i) = tmp(1:3);
    contr(:,i) = res{i,3};
    norm_lambda(i) = norm(res{i,2});
end

norm_t = zeros(1,n_timepoints);
for i=1:n_timepoints
norm_t(i) = norm(cam_pos(:,i) - pos(:,i));
end
save('visualization/ExampleSkier/skierCamPos_Short.mat');

%% Plot functions
% Plot cost function
load('skierCamPos_Short.mat');

relpath = 'visualization/ExampleSkier/';

figure;
plot(costF,'r')
title('Cost function value');
print([relpath, 'costF'], '-dsvg');

% Plot CamPos and Pos
figure;
plot3(cam_pos(1,:), cam_pos(2,:), cam_pos(3,:),'b')
hold on
plot3(pos(1,:), pos(2,:), pos(3,:),'r')
view(-116,16);
title('Flight in R^3');
legend({'Given camera position' , 'Drone'}, 'Location', 'northwest');
print([relpath, 'r3Plot'], '-dsvg');

%Plot
figure
plot(norm_t,'r');
title('Distance between camPosition and computed position');
print([relpath, 'norm_t'], '-dsvg');

% Plot controls
figure
subplot(2,1,1);
title('Controls');
hold on
for k = 1:2
    contrPlot = subplot(2,1,k);
    hold on
    plot(contr(k,:)', 'b');
    plot(contr(k+2,:)', 'r');
    legend(contrPlot, {['Engine ' , int2str(k)] ,['Engine ' , int2str(k+2)]}, 'Location', 'north','FontSize', 12);
end
print([relpath, 'controls'], '-dsvg');

%Plot norm lambdas
figure
plot(norm_lambda', 'r');
title('Norm of lambdas');
print([relpath, 'norm_lambda'], '-dsvg');

%% For Screencast of downhill setting
close all;
pause(1);

hLarge = figure('name', 'Quadrocopter following a skier');
% Set size of window to fullHD
set(hLarge , 'Position', [0 0, 1920 , 1080 ]);


pause(4);

for i = 1:n_timepoints
    r3plot = subplot(2,3,[1,2,4,5]);
   
    plot3(cam_pos(1,1:i), cam_pos(2,1:i), cam_pos(3,1:i),'b')
    hold on
    plot3(pos(1,1:i), pos(2,1:i), pos(3,1:i),'r');
    plot3(cam_pos(1,i), cam_pos(2,i), cam_pos(3,i),'b.', 'markersize', 19);
    plot3(pos(1,i), pos(2,i), pos(3,i),'r.', 'markersize', 24);
    axis([-20 20  -5 305 -5 275]);
    legend(r3plot, {'Given camera position' , 'Drone'}, 'Location', 'northwest', 'FontSize', 14);
    set(gca, 'FontSize', 12);
    view(-116,16);
    title( 'Flight in R^3','FontSize', 14);
    hold off
    
    subplot(2,3,3);
    plot(norm_t(1:i),'r');
    hold on 
    plot( i, norm_t(i), 'r.', 'markersize', 30);
    axis( [0 n_timepoints 0 15] );
    set(gca, 'FontSize', 12);
    title( 'Distance between actual position and given camera position','FontSize', 12);
    hold off
    
    
    subplot(2,3,6);
    plot(costF(1:i),'r');
    hold on 
    plot(i,costF(i), 'r.', 'markersize', 30);
    axis( [0 n_timepoints 0 220])
    set(gca, 'FontSize', 12);
    title( 'Cost function','FontSize', 12);
    hold off
    
    
    pause(0.01);
end

print([relpath, 'Animation'], '-dsvg');

%% Show animation in nice framework

posi = zeros(n_timepoints,3);
ang = zeros(n_timepoints,4);
dot_pos = zeros(n_timepoints,3);
dot_ang = zeros(n_timepoints,3);

for i = 1:n_timepoints
   tmp  = res{i,1};
    posi(i,:) = tmp(1:3)';
    ang(i,:) = tmp(4:7)';
    dot_pos(i,:) = tmp(8:10)';
    dot_ang(i,:) = tmp(11:13)';
end

angle = zeros(n_timepoints, 3);
[angle(:, 1), angle(:, 2), angle(:, 3)] = quat2angle((ang));
% Prepare variables
global visualization

visualization.yout = [ dot_ang , angle , dot_pos , posi, rand(n_timepoints, 3), norm_t'];
visualization.tout = ones(n_timepoints,1);

QuadAnim41;