classdef Constraints < GenConstraints & TestEnv
    % classConstraints providing discretized ODE constraint using forward euler
    
    properties %Realtime Properties
        activeSet;
        n_addConstr = 8;
    end
    properties %
        flagWind;
        flagWindD;
        Wind;
        WindD;
    end
    methods
        %constructor
        function cC = Constraints(varargin)
            
            dode = [];
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
            elseif (nargin == 1)
                if (isa(varargin{1},'MultiShooting'))
                    dode = varargin{1};
                else
                    error('wrong class type for discretized ode');
                end
            else
                error('wrong number of inputs');
            end
            cC@GenConstraints(dode);
            cC.flagWind = true;
            cC.flagWindD  = true;
        end
        
        
        %Check
        
%         % other functions
%         function [ineq_con,eq_con,ineq_conD,eq_conD] = constr(obj)
%             % provide the equality and inequality constraints (along
%             % with their jacobians by calling get_eq_con and
%             % get_eq_conD
%             ineq_con                = obj.get_ineq_con();
%             [eq_con, eq_conD]       = obj.get_eq_con();
%             ineq_conD               = obj.get_ineq_conD();
%         end %Check
        % InEqualities
        function ineq_con = get_ineq_con(obj)
            ineq_con = [];%obj.InEqCon;
        end %Check
        
        function ineq_conD = get_ineq_conD(obj)
            ineq_conD = [];%obj.InEqConD';
        end %Check
        
        function ineq_conDD = get_ineq_conDD(obj)
            ineq_conDD = [];%obj.InEqConDD;
        end
        
        % Equalties
        function [eq_con] = get_eq_con(obj) 
            % the equality constraint of the ocp
            % combine the discretized ode with the boundary conditions
            
            %xbc         = obj.dyn.environment.xbc;
            state_mat   = obj.dyn.state;
            
            
            [H] = obj.dode.h();
            
            
            
            eq_con      = [H; 
                           obj.getWind(state_mat);
                          ];
        end 
        
        function eq_conD = get_eq_conD(obj)
            % the Jacobian of the equality contraints of the ocp
            [r,q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
            n_var = 17;
            
            
            
            
           [H, hD] = obj.dode.h();
            
            eq_conD     = [hD;  ...
                obj.spWindD();
                %obj.EqConD ...
                ]';
            
        end
        
        function eq_conDD = get_eq_conDD(varargin)
            % the Hessian of the equality constraints of the ocp.
            % Returns a cell array of Hessians if only the object
            % is provided. If the object and the Langrange multipliers are
            % provided, the function returns the Hessian of the Lagrangian
            if (nargin == 1)
                obj         = varargin{1};
                
                [r, q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
                n_timepoints = obj.dyn.environment.n_timepoints;
                
                eq_conDD = [obj.dode.hDD(); cell((n_timepoints-1)*n_state,1)];
                for i=0:(n_timepoints-1)*n_state-1
                    eq_conDD{end-i} = sparse((n_int+1)*(n_state+n_contr),...
                        (n_int+1)*(n_state+n_contr));
                end
                eq_conDD = [eq_conDD]; %obj.EqConDD];
            elseif (nargin == 2)
                lambda      = varargin{2};
                H           = varargin{1}.dode.hDD();
                
                if (~isempty(lambda))
                    eq_conDD    = {lambda(1)*H{1}};
                    for i=2:length(lambda)
                        eq_conDD{1} = eq_conDD{1} + lambda(i)*H{i};
                    end
                else
                    eq_conDD = cell(1,1);
                end
            else
                eq_conDD = cell(1,1);
                %             else
                %                 error('Error: Not yet implemented');
            end
        end
        
        function H = get_eqhess(varargin)
            % interfacing get_eq_conDD
            if (nargin == 1)
                H = get_eq_conDD(varargin{1});
            elseif (nargin == 2)
                H = get_eq_conDD(varargin{1},varargin{2});
            else
                H = [];
            end
        end
        
        function H = get_ineqhess(varargin)
            % interfacing get_eq_conDD
            if (nargin == 1)
                H = get_ineq_conDD(varargin{1});
            else
                H = [];
            end
        end
        
    end
    methods %Realtime
        function eq_con = get_eq_con_at_t(obj,t)
            % GET_EQ_CON_AT_T Returns all equality constraints for the realtime
            % optimization problem at time t.
            
            state_mat   = obj.dyn.state;
            s_t = state_mat(:,1);
            multShoot = obj.dode.h();
            
            eq_con      = [obj.dyn.environment.wind(s_t, t) - state_mat(:,1);
                multShoot();
                ];
        end
        
        %TODO: get_eq_con_at_t und get_ineq_con_at_t ist inkonsistent.
        function ineq_con = get_ineq_con_at_t(obj,t)
            % GET_INEQ_CON_AT_T Returns all inequality constraints at
            % timepoint t.
            ineq_con =  obj.ineqCon();
            ineq_con = ineq_con( (t-1) * obj.n_addConstr +1 : t * obj.n_addConstr );
        end
        
        function eq_conD = get_eq_conD_block_t(obj,t)
            % the Jacobian of the equality contraints of the ocp
            [r,q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
            [h, multShootD] = obj.dode.h();
            
            eq_conD     = multShootD( (t-1) * n_state +1  : t*n_state, (t-1)* (n_contr+ n_state) +1 : t  * (n_contr+ n_state));
            
            %eq_conD = [ A,B];
        end
        
        function ineq_conD = get_ineq_conD_block_t(obj,t)
            ineq_conD = obj.ineqConD_at_t(t);
            
            %ineq_conD = [ 0, D];
        end
        
        function eq_conDD = get_eq_conDD_at_t(o, lambda, t)
            % GET_EQ_CONDD_AT_T Returns the second derivative of the eq_constraints for timepoint t
            n_state = length(lambda);
            if isa(o.dode.solver, 'ForwEuler')
                
                H = o.dode.hDD();
                
                if (~isempty(lambda))
                    eq_conDD = lambda(1) * H{ (t-1) * n_state +1 };
                    for i = 2:n_state
                        eq_conDD = eq_conDD + lambda(i) * H{ (t-1) * n_state + i};
                    end
                else
                    error('lambda should not be empty');
                end
                
            else
                % derive numerically, as implementation is too time
                % consuming
                H = cell(n_state, 1);
                for i = 1:n_state
                    func = @() o.numDerivHelper(i,t);
                    H{i} = o.numDiff_nxnD(t,func);
                    disp(['State No. :', int2str(i)]);
                end
                
                eq_conDD = 0;
                for  i = 1:n_state
                    eq_conDD = eq_conDD + lambda(i) * reshape(H{i}, size(H{i},2),size(H{i},3));
                end
                
            end
        end
        
        function ineqCon = ineqCon(o)
            % INEQCON Calculates the inequality constraints, such that
            % for every control signal q_i holds: u_min <= q_i <= u_max
            n_intervals = o.dyn.environment.n_intervals;
            ineqCon = zeros(o.n_addConstr*(n_intervals+1),1);
            
            for i = 1:n_intervals+1
                ineqCon( (i-1) * o.n_addConstr +1: i * o.n_addConstr) = ...
                    [ o.dyn.vec( (i-1) * o.dyn.robot.n_var + o.dyn.robot.n_state +1 : i * o.dyn.robot.n_var) - o.dyn.robot.u_max; ...
                    o.dyn.robot.u_min - o.dyn.vec( (i-1) * o.dyn.robot.n_var + o.dyn.robot.n_state +1 : i * o.dyn.robot.n_var)];
            end
        end
        
        function ineqConD = ineqConD_at_t(o,t)
            % INEQCOND_AT_T Calculats the derivative of ineqCon for
            % timestep t. Which ensures, that every control q_i is in a
            % given interval. As the result is not timedependent, t is not
            % used.
            ineqConD = [sparse(o.dyn.robot.n_contr, o.dyn.robot.n_state) , speye(o.dyn.robot.n_contr)];
            ineqConD = [ineqConD; -ineqConD];
        end
        
        function mu = checkIfActive(o,mu)
            % CHECKIFACTIVE Checks which constraint is active and which isnt and stores it
            % in the activeSet property.
            
            bIneq = o.ineqCon() >= 0;
            bMu = mu > 0;
            %Check mus from previous iteration
            o.activeSet = bIneq & bMu;
            %Set inactive mus to 1
            mu(~o.activeSet) = 1;
        end
        
        function activeSet_k = getActiveSet(o,k)
            % GETACTIVESET Gives you a boolean vector of size o.n_addConstr
            % indicating, which additional inequality constraint is active
            % and which isn't.
            activeSet_k = o.activeSet( (k-1) * o.n_addConstr +1 : k * o.n_addConstr);
        end
        
        function hD =  numDerivHelper(obj,i,t)
            n_state = obj.dode.dyn.robot.n_state;
            
            [h , hDhelp ] = obj.dode.h();
            hD = hDhelp( (t-1) * n_state + i, : )';
        end
        
        
        
        
    end
    methods %Helper
        
        function res = getWind(obj, states)
            
            n_state = 13;
            n_timepoints = obj.dyn.environment.n_timepoints;
            if obj.flagWind %Wind nur einmal ueberall Zeitraueme initialisieren
                 
                obj.Wind = zeros(n_state * (n_timepoints-1), 1);

                for timepoint = 1:(n_timepoints-1)
                    st_ = zeros(n_state, 1);
                    obj.Wind( (timepoint-1) * n_state + 1: timepoint * n_state )= obj.dyn.environment.wind(st_, timepoint);
                end
                obj.flagWind = false;
            end
            res = zeros(n_state * (n_timepoints-1), 1);
            for timepoint = 1:(n_timepoints-1)
                res((timepoint-1) * n_state + 1: timepoint * n_state ) = states(:, timepoint ) + obj.Wind( (timepoint-1) * n_state + 1: timepoint * n_state );
            end
            
        end
        function res = spWindD(obj)
            if obj.flagWindD
                
                [r,q,v,omega,u,Iges,IM,m,kT,kQ,d,g, n_int, n_state, n_contr] = getParams(obj);
                n_timepoints = obj.dyn.environment.n_timepoints; 
                n_var = 17;
                sp = zeros((n_timepoints-1)*n_state, 3);

                for timepoint = 1:(n_timepoints-1)
                    sp((timepoint-1)*n_state +1:timepoint*n_state, 1) = (timepoint-1)*n_state +1:timepoint*n_state;
                    sp((timepoint-1)*n_state +1:timepoint*n_state, 2) = (timepoint-1)*n_var +1:(timepoint-1) * n_var + n_state;
                    sp((timepoint-1)*n_state +1:timepoint*n_state, 3) = ones(1, n_state);
                end

                obj.WindD = sparse(sp(:,1),sp(:,2),sp(:,3),...
                    (n_timepoints-1)*n_state,(n_timepoints)*(n_var));
                obj.flagWindD = false;
            end
            res = obj.WindD;
        end
    end
    methods %TestHelper
        
        %Some help functions (typically overwritten in subclasses)
        function [vec_old, n,m, n_timepoints, dyn] = setup(obj,func)
            vec_old = obj.dyn.vec;
            n_timepoints = obj.dyn.environment.n_timepoints;
            dyn = obj.dyn;
            n = obj.dyn.robot.n_var;
            m = size(func());
            m=m(1);
        end
        
    end    
    methods(Test)
%Check
        function test_get_eqjac(obj)
            %TEST_GET_EQJAC This method derives numerically get_eqfunc and compares it
            %with get_eqjac

            n_intervals = uint16(10);
            obj.setupTest(n_intervals, ForwEuler());

            func = @() obj.get_eq_con;
            
            numDiff = obj.numDiff_nD_AllT(func);
            
            anaDiff = obj.get_eq_conD();
            anaDiff = anaDiff';
            
            obj.assertSize(anaDiff, size(numDiff) );
            obj.assertLessThan(max(abs(anaDiff - numDiff)), obj.tol);

        end
        
% Check
        function test_get_eqhess(obj)
            % TEST_GET_EQHESS This method derives numerically get_eqjac and
            % compares it with get_eqhess
            n_intervals = uint16(2);
            obj.setupTest(n_intervals,ForwEuler());

            func = @() obj.get_eq_conD()'; %TODO: passt das?
            anaDiff = obj.get_eqhess();
            numDiff = obj.numDiff_nxnD_AllT(func);

            size_nDiff_i = (obj.dyn.robot.n_var) * (n_intervals +1 );
            for i = 1:length(anaDiff)
                numDiff_i = reshape(numDiff(i,:,:), [size_nDiff_i size_nDiff_i]);
                obj.assertSize(anaDiff{1}, size(numDiff_i));
                obj.assertLessThan(max(abs(anaDiff{i} - numDiff_i)), obj.tol);
            end
         end

    end
    
    methods
        function setupTest(obj, n_intervals,integrator)
            
            n_int_ = n_intervals;
            % Quadrocopter soll 5 Meter hoch fliegen
            xbc = [         ... Variablenname L�nge   Name
                ... Anfangsbedingung
                0, 0, 0.5,  ...     r           3      Ortsvektor
                1, 0, 0, 0, ...     q           4      Quaternion (Einheitsquaternion)
                0, 0, 0,    ...     v           3      Translatorische Geschwindigkeit
                0, 0, 0;    ...     w           3      Winkelgeschwindigkeit
                ... Endbedingung
                0, 0, 0,    ...
                1, 0, 0, 0, ...
                0, 0, 0,    ...
                0, 0, 0     ...
                ];
            
            env = Environment();
            env.xbc = xbc;
            env.wind = @(s_t ,t ) s_t + 0.1 * [rand(3,1); zeros(10,1)];
            env.setUniformMesh(n_int_);
            
            robot = Quadrocopter();
            

            dyn_ = BasisQDyn(robot, env, integrator);
            dyn_.vec = rand(17 * (n_int_ + 1), 1);
            
            obj.dode = MultiShooting(dyn_);
            obj.dyn = obj.dode.dyn;
            obj.dode.noCaching = true;
            
        end
        
        
    end
end