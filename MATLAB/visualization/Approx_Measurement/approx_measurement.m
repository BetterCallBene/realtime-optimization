clear all


horizon = 15;

cLagrange = Lagrange();
app1 = @(cRTSolver,t) cLagrange.getLDD_approx_costDDpAlphaI(cRTSolver, t, zeros(17,1) ) ;
app2 = @(cRTSolver,t) cLagrange.getLDD_approx_costDDpAlphaI(cRTSolver, t, 0.01 * ones(17,1) ) ;
app3 = @(cRTSolver,t) cLagrange.getLDD_approx_costDDpAlphaI(cRTSolver, t, 0.1 * ones(17,1) ) ;
app4 = @(cRTSolver,t) cLagrange.getLDD_approx(cRTSolver,t);
app5 = @(cRTSolver,t) cLagrange.getLDD(cRTSolver,t);


LDDcell = { app5 ,app4, app3, app2, app1};



n_rep = 10;

time = zeros(length(LDDcell), n_rep);
costF = cell(length(LDDcell), n_rep);
error = cell(length(LDDcell), n_rep);
n_timepoints = 20 ; %How many timepoints, do we want to calculate.

pos = zeros(3,n_timepoints);
global TEST
TEST = true;
cCost = CostsComplet();
TEST = false;
for i = 1:n_timepoints
    pos(:,i) = [0;0;0];
end

rand_hF = rand(3,1);
rand_hM = rand(3,1);



%% Get time values 
for h = 1:length(LDDcell)
    for j = 1:n_rep
        
        %%Setup
        
        pointPerSecond = 1;
        
        env = Environment();
        env.horizon = horizon;
        %Die Dynamik wird nur auf dem Horizon betrachtet:
        n_intervals = env.setUniformMesh1(horizon+1,pointPerSecond);
        
        cQ = Quadrocopter();
        
        % Wahl des Integrators
        tol = 1e-2;
        opts = odeset('RelTol',tol,'AbsTol',0.1*tol);
%         cIntegratorExt = ode15sM(opts);
        cIntegratorExt = ForwEuler(opts);
        
        cQExt = QuadrocopterExt(cQ, env, cIntegratorExt);
        cQExt.steadyPoint = [];  %steadyPoint initialisieren: SteadyPoint ist eine globale Variable!!
        %cQExt.hForceExt = @(v) 0.5 * rand_hF + cQ.getF_w(v);
        %cQExt.hMomentExt = @() 0.5 * rand_hM;
        %Neue Windfunktion
        %env.wind = @(s_t, ctr)  cQExt.wind(s_t, ctr);
         env.wind = @(s_t ,t ) s_t + 0.01 * [rand_hF; zeros(10,1)];
%        env.wind = @(s_t ,t ) s_t ;
        % Initialisierung der Dynamik
        cBQD = BasisQDyn(cQ, env, cIntegratorExt);
        
        % Initialisierung des Multiple Shootings
        cMultShoot = MultiShooting(cBQD);
        
        % Initialisierung der Nebenbedingungen
        cConst = Constraints(cMultShoot);
        
        % Initialisierung Kostenfunktion
        cCost = CostsComplet(cBQD, 1, 0.5, 1, 1);
        
        %Define Cam Position function
        cCost.cam_pos = @(t) [0;0;0];
        
       % Choose starting values
        
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
        getLDD = LDDcell{h};
        
        t_tic = tic;
        %Use realtime solver
        [res, est_y  ] = cRTSolver.fminrt(getLD, getLDD, n_timepoints);
        time(h,j)  =  toc(t_tic);     
        
         tmp = zeros(n_timepoints,1);
        tmp2 = zeros(n_timepoints,1);
        for i = 1:n_timepoints
            mtp = res{i,1};
            tmp2(i) = norm(mtp(1:3) - pos(:,i));
            tmp(i) = res{i,5};
        end
        costF{h,j} = tmp;
        error{h,j} = tmp2;
    end
end


%% Eval

errMin = -ones(length(LDDcell),1);
coFMin = -ones(length(LDDcell),1);

timeMin = min(time');

for h = 1:length(LDDcell)
    for l = 1:n_rep
         actError = error{h,l};
         actCost = costF{h,l}; 
        if( errMin(h) < 0  ||  errMin(h) > actError(end)  )
            errMin(h)  = actError(end);
        end
        
        if ( coFMin(h) < 0 || coFMin(h) > actCost(end) )
           coFMin(h) = actCost(end); 
        end
        
    end
    
end


%plot(timeMin);
% 
% figure;
% 
% plot(errMin);
% figure;
% 
% plot(coFMin);
% 

figure
val = [timeMin' , errMin , coFMin];


val = val([1,3,4,5], :);
names = {'Exact', 'costDD+0.1', 'costDD+0.01', 'costDD'};
bar(val);



subplot(1,3,1)
bar(val(:,1)-6.5,'red')
title('Time');
set(gca,'XTickLabel',names, 'XTick', 1:4); 
% set(gca, 'FontSize', 10);


subplot(1,3,2)
bar(val(:,2)-0.027,'blue')
title('Distance');
set(gca,'XTickLabel',names, 'XTick', 1:4); 
% set(gca, 'FontSize', 10);


subplot(1,3,3)
bar(val(:,3)-13,'red')
title('Cost function');
set(gca,'XTickLabel',names, 'XTick', 1:4); 
% set(gca, 'FontSize', 10);


relpath = 'visualization/Approx_Measurement/';
print([relpath, 'approxPlot'], '-dsvg');

