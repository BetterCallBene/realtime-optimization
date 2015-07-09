clear all

horizons = [4, 15, 40];
n_rep = 7;

time = zeros(length(horizons), n_rep);
costF = cell(length(horizons), n_rep);
error = cell(length(horizons), n_rep);
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
for h = 1:length(horizons)
    for j = 1:n_rep
        
        %%Setup
        horizon = horizons(h);
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
        %cQExt.hForceExt = @(v) 0.5 * rand_hF + cQ.getF_w(v);
        %cQExt.hMomentExt = @() 0.5 * rand_hM;
        %Neue Windfunktion
        %env.wind = @(t, s_t, ctr)  cQExt.wind(s_t, ctr);
        env.wind = @(t, s_t, ctr) s_t + 0.5 * [rand_hF; zeros(10,1)];
        % Initialisierung der Dynamik
        cBQD = BasisQDyn(cQ, env, cIntegrator);
        
        % Initialisierung des Multiple Shootings
        cMultShoot = MultiShooting(cBQD);
        
        % Initialisierung der Nebenbedingungen
        cConst = Constraints(cMultShoot);
        
        % Initialisierung Kostenfunktion
        
        cCost = CostsComplet(cBQD, 1, 0.1, 1, 1);
        
        
        
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
        getLDD = @(cRTSolver,t) cLagrange.getLDD_approx_costDDpAlphaI(cRTSolver, t, zeros(17,1) ) ;
        
        
        t_tic = tic;
        %Use realtime solver
        [res, est_y  ] = cRTSolver.fminrt(getLD, getLDD, n_timepoints);
        time(h,j)  =  toc(t_tic);       
        
    end
end



%% Get error values with wind and aerodynamics

for h = 1:length(horizons)
    for j = 1:n_rep
        
        %Setup
        horizon = horizons(h);
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
        cQExt.hForceExt = @(v) 0.5 * rand_hF + cQ.getF_w(v);
        cQExt.hMomentExt = @() 0.5 * rand_hM;
        %Neue Windfunktion
        env.wind = @(t, s_t, ctr)  cQExt.wind(s_t, ctr);
        %env.wind = @(t, s_t ,ctr ) s_t + 0.5 * [rand(3,1); zeros(10,1)];
        % Initialisierung der Dynamik
        cBQD = BasisQDyn(cQ, env, cIntegrator);
        
        % Initialisierung des Multiple Shootings
        cMultShoot = MultiShooting(cBQD);
        
        % Initialisierung der Nebenbedingungen
        cConst = Constraints(cMultShoot);
        
        % Initialisierung Kostenfunktion
        
        cCost = CostsComplet(cBQD, 1, 0.1, 1, 1);
        
        
        
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
        getLDD = @(cRTSolver,t) cLagrange.getLDD_approx_costDDpAlphaI(cRTSolver, t, zeros(17,1) ) ;
        
        % Calculate the solution with fminrt
               
        %Use realtime solver
        [res, est_y  ] = cRTSolver.fminrt(getLD, getLDD, n_timepoints);
        
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

%% Evaluation


subplot(length(horizons),n_rep,1)

for k = 0:length(horizons)-1
    for l = 1:n_rep
    subplot(length(horizons),n_rep,k*n_rep + l );
 %   ax2 = axes('YAxisLocation','right');
    plot(error{k+1,l});
    title(['Horizon: ', int2str(horizons(k+1)), ' Time: ', num2str(time(k+1,l)/n_timepoints)]);
    axis([ 1 n_timepoints 0 (max(error{1,1})+1) ] );
    %ax2 = axes('YAxisLocation','right');
   % set(gca, 'ytick', (linspace(0, max(error{1,1})+1 , 5)));
    end
end

figure;

endError = zeros(length(horizons),1);
endT = endError;
names = { };
for k = 0:length(horizons)-1
    for l = 1:n_rep
        actError  = error{k+1,l};
        if( endError(k+1) < actError(end))
            endError(k+1) = actError(end);
        end
    end
    endT(k+1) = min(time(k+1,:)) / n_timepoints;
    names = { names{1:end} , ['Horizon: ' , int2str(horizons(k+1))] };
end

figure;
ende = [ endT, endError];
b = bar(ende);
set(gca,'XTickLabel',names, 'XTick', 1:3 ); 

b(1).FaceColor = 'red';
b(2).FaceColor = 'blue';
             
relpath = 'visualization/Horizon_Measurement/';
print([relpath, 'horizonPlot'], '-dsvg');
legend('Time in s', 'Error in m');
set(gca, 'FontSize', 12);
