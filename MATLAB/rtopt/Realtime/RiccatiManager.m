classdef RiccatiManager < TestEnv
    % RICCATIMANAGER Central class to perform a riccati solution
    
    properties(Access=private)
        robot;
        A;
        B;
        M;
        P;
        Q;
        R;
        
        D;
        
        nabla_s_star;
        nabla_lambda;
        nabla_q;
        
        nabla_mu;
        
        z3;
        z4;
        z5;
        
        horizon;
        n_lambda;
        n_state;
        n_contr;
        n_mu; % Das ist eine Cell, da n_mu die aktiven Constraints in jedem Zeitschritt unterschiedlich sind.
        n_var; % Das ist eine Cell, da n_mu die aktiven Constraints in jedem Zeitschritt unterschiedlich sind.
    end
    
    properties
        delta_lambda;  
        delta_s;
        delta_q;
        delta_mu;
    end
    
    
    methods
        function o = RiccatiManager(varargin)
            
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                else
                    o.horizon = 20;
                    o.robot = Quadrocopter();
                end
            elseif(nargin == 1)
                o.horizon = varargin{1};
                o.robot = Quadrocopter();
            elseif(nargin == 2)
                o.horizon = varargin{1};
                if (isa(varargin{2},'Model'))
                    o.robot = varargin{2};
                end
            else
                error('wrong number of inputs');
            end
            
            o.initialize();
        end
        
        function doStep(o,i, LDD_i, LD_i, n_mu_i )
            o.n_mu{i} = n_mu_i;
            o.Q{i} = LDD_i(1:o.n_state,1:o.n_state);
            o.nabla_lambda{i} = LD_i(1:o.n_lambda);
            o.n_var{i} = o.n_lambda + o.n_state + o.n_contr + n_mu_i;
            
            if(i == o.horizon+1)
                %Dann ist n_mu_i = 0 da dann keine controls vorhanden sind
                o.P{i} = o.Q{i};
                o.nabla_s_star{i} = LD_i(o.n_lambda+1:o.n_lambda + o.n_state);
                
            else
                o.M{i} = LDD_i(1:o.n_state,o.n_state +1 :o.n_state + o.n_contr);
                o.R{i} = LDD_i(o.n_state +1 :o.n_state + o.n_contr, o.n_state+1 :o.n_state + o.n_contr);
                o.A{i} = LDD_i(o.n_state + o.n_contr + n_mu_i + 1: o.n_var{i}, 1:o.n_state);
                o.B{i} = LDD_i(o.n_state + o.n_contr + n_mu_i + 1: o.n_var{i}, o.n_state +1 : o.n_state + o.n_contr);
                o.nabla_q{i} = LD_i(o.n_lambda + o.n_state +1  : o.n_lambda + o.n_state + o.n_contr);
                if(n_mu_i == 0)
                    % Benutze Regel aus Riccati_woConstr
                    
                    o.P{i} = o.Q{i} + (o.A{i}' * o.P{i+1} * o.A{i}) - ...
                        (o.M{i} + o.A{i}' * o.P{i+1} * o.B{i}) * ...
                        ((o.R{i} + o.B{i}' * o.P{i+1}  * o.B{i}) \ ...
                        (o.M{i}'  + o.B{i}' * o.P{i+1} * o.A{i}));
                    
                    o.nabla_s_star{i} = LD_i(o.n_lambda +1 : o.n_lambda + o.n_state) + ...
                        o.A{i}' * o.P{i+1} * o.nabla_lambda{i+1} + ...
                        o.A{i}' * o.nabla_s_star{i+1} - ...
                        (o.M{i} + o.A{i}' * o.P{i+1} * o.B{i}) * ...
                        (( o.R{i} +  o.B{i}' * o.P{i+1} *   o.B{i}) ...
                        \ ( o.nabla_q{i} + o.B{i}' * o.P{i+1} *  o.nabla_lambda{i+1} + o.B{i}' * o.nabla_s_star{i+1}));
                    %TODO: hier kann man o.B{i}' und auch o.A{i}' ausklammern, machen wir
                    %aber erst, wenn der Test geht.
                else
                    % Benutze neue Regel mit Constraints
                    o.D{i} = LDD_i(o.n_state + o.n_contr + 1 : o.n_state + o.n_contr + n_mu_i ,  o.n_state +1 : o.n_state + o.n_contr);
                    o.nabla_mu{i} = LD_i( o.n_lambda + o.n_state + o.n_contr +1 : o.n_lambda + o.n_state + o.n_contr + n_mu_i);
                    
                    o.z4{i} = [o.M{i}'  + o.B{i}' * o.P{i+1} * o.A{i}; zeros(n_mu_i, o.n_state)];
                    o.z3{i} = [ o.R{i} +  o.B{i}' * o.P{i+1} *   o.B{i}, o.D{i}' ;...
                        o.D{i} , zeros(n_mu_i)];
                    z2 = [o.M{i} + o.A{i}' * o.P{i+1} * o.B{i}, zeros(o.n_state, n_mu_i)];
                    z1 =  o.Q{i} + (o.A{i}' * o.P{i+1} * o.A{i});
                    
                    o.P{i} = z1 - z2 * (o.z3{i} \ o.z4{i});
                    
                    o.z5{i} = [ o.nabla_q{i} + o.B{i}' * ( o.P{i+1} * o.nabla_lambda{i+1} + o.nabla_s_star{i+1}) ; ...
                        o.nabla_mu{i}];
                    
                    z6 = LD_i(o.n_lambda +1 : o.n_lambda + o.n_state) + ...
                        o.A{i}' * (o.P{i+1} * o.nabla_lambda{i+1} + o.nabla_s_star{i+1});
                    
                    o.nabla_s_star{i} = z6 - z2 * (o.z3{i} \ o.z5{i});
                    
                end
            end
        end
        
        function solveStep(o,i)
            if ( i ==1 )
                o.delta_lambda{i} = -o.P{i}* o.nabla_lambda{i} - o.nabla_s_star{i};
                o.delta_s{i} = -o.nabla_lambda{i};
            else
                z1 = (o.A{i-1} * o.delta_s{i-1} + o.B{i-1} * o.delta_q{i-1});
                
                % solve delta lambda
                o.delta_lambda{i} = -o.P{i}* o.nabla_lambda{i} - o.nabla_s_star{i} + o.P{i} * z1 ;
                % solve delta s
                o.delta_s{i} = - o.nabla_lambda{i}  + z1;
            end
            
            %solve for delta_q and delta_mu
            if( i ~= o.horizon +1)
                if ( o.n_mu{i} == 0)
                    o.delta_q{i} = (o.R{i} + o.B{i}' * o.P{i+1} * o.B{i}) \ ...
                        ( ...
                        o.nabla_q{i} + ...
                        o.B{i}' * o.P{i+1} * o.nabla_lambda{i+1} + ...
                        o.B{i}' * o.nabla_s_star{i+1} - ...
                        ( o.M{i}' + o.B{i}' * o.P{i+1} * o.A{i} ) * o.delta_s{i} );
                    
                else
                    z6 = o.z3{i} \ (o.z5{i} - (o.z4{i} * o.delta_s{i}));
                    
                    o.delta_q{i} = z6(1:o.n_contr);
                    o.delta_mu{i} = z6(o.n_contr +1 : o.n_contr + o.n_mu{i});
                    
                end
            end
        end
        
        function delta_mu = assembleMu(o,activeSet_i,i)
            % ASSEMBLEMU builds a mu, which also considers the inactive constraints
            delta_mu = zeros(length(activeSet_i),1);
            delta_mu(activeSet_i) = o.delta_mu{i};
        end
    end
    
    methods(Test)
        
        function testRiccati_0(o)
            o.doTestWithHorizon(0,0,1e-8);
        end
        
        function testRiccati_1(o)
            o.doTestWithHorizon(1,0,1e-8);
        end
        
        function testRiccati_20(o)
            o.doTestWithHorizon(20,0,1e-8);
        end
        
        function testRiccati_100(o)
            o.doTestWithHorizon(100,0,1e-8);
        end
        
        function testRiccatiwConstr_0(o)
            o.doTestWithHorizon(0,randi(4,1),1e-9);
        end
        
        function testRiccatiwConstr_1(o)
            o.doTestWithHorizon(1,randi(3,1),1e-7);
        end
        
        function testRiccatiwConstr_20(o)
            o.doTestWithHorizon(20,2,1e-4);
        end
        
    end
    
    methods
        
        function initialize(o)
            %Initialize storage
            o.n_lambda = o.robot.n_state;
            o.n_state = o.robot.n_state;
            o.n_contr = o.robot.n_contr;
            o.n_mu = cell(o.horizon+1,1);
            o.n_var = cell(o.horizon+1,1);
            
            o.A = cell(o.horizon,1);
            
            
            o.B = cell(o.horizon,1);
            o.M = cell(o.horizon,1);
            o.P = cell(o.horizon+1,1);
            o.Q = cell(o.horizon+1,1);
            o.R = cell(o.horizon,1);
            o.D = cell(o.horizon,1);
            
            o.z3 = cell(o.horizon,1);
            o.z4 = cell(o.horizon,1);
            o.z5 = cell(o.horizon,1);
            
            o.nabla_s_star = cell(o.horizon+1, 1);
            o.nabla_lambda = cell(o.horizon+1,1);
            o.nabla_q = cell(o.horizon ,1);
            o.nabla_mu = cell(o.horizon, 1);
            
            %    o.delta = zeros( o.n_var * (o.horizon +1 ) - (o.n_mu +  o.n_contr), 1);
            
            
            o.delta_lambda = cell(o.horizon +1, 1);
            o.delta_s = cell(o.horizon +1, 1);
            o.delta_q = cell(o.horizon, 1);
            o.delta_mu = cell(o.horizon,1);
        end
        
        function doTestWithHorizon(o, hori,n_addConstr, tolerance)
            o.horizon = hori;
            o.robot = Quadrocopter();
            o.initialize();
            
            LDDi = zeros(30 + n_addConstr);
            %Find an invertible matrix
            eigLDDi = abs(eig(LDDi));
            while (min(eigLDDi) < 1e-4 || max(eigLDDi) > 1e4)
                LDDi = 10 * eye(30+n_addConstr) +  2*rand(30+n_addConstr);
                LDDi(18:end, 18:end) = zeros(13+n_addConstr);
                
                %Set derivatives of addConstr after states to 0
                LDDi(17+1: 17+n_addConstr, 1:13) = zeros(n_addConstr, 13);
                LDDi(13+4+1:13+4+n_addConstr, 1:13) = zeros(n_addConstr,13);
                LDDi(1:13,13+4+1:13+4+n_addConstr) = zeros(n_addConstr,13)';
              
                %LDDi must be symmetric
                LDDi = .5* (LDDi + LDDi');
                eigLDDi = abs(eig(LDDi));
            end
            getLDD = @(i) LDDi;
            hesse_L = o.buildUpHesse(getLDD,30+n_addConstr);
            deltay = 10 * rand((o.horizon+1) * (30 + n_addConstr) - (4+n_addConstr) ,1);
            grad_L = hesse_L * deltay;
            
            disp(' ');
            tic;
            %Build up Riccati stack
            for i = o.horizon+1:-1:1
                if i == o.horizon +1
                    o.doStep(i,LDDi, grad_L( (i-1) * (30+n_addConstr) +1 : end ,1), n_addConstr );
                else
                    o.doStep(i,LDDi, grad_L( (i-1) * (30+n_addConstr) +1 : i * (30+n_addConstr),1), n_addConstr);
                end
            end
            
            % Solve the the stack
            for i = 1:o.horizon +1
                o.solveStep(i);
            end
            
            timeRiccati = toc;
            disp(strcat('Time for Riccati with Horizon',{' '}, int2str(hori), ':'));
            disp(num2str(timeRiccati));
            
            %Compare calculation time with \
            tic;
            hesse_L = o.buildUpHesse(getLDD,30+n_addConstr);
            delta_matlab = hesse_L \ grad_L;
            timeMatlab = toc;
            disp(strcat('Time for \ with Horizon',{' '}, int2str(hori), ': '));
            disp(num2str(timeMatlab));
            
            
            %Check if result is correct
            delta = o.assembleDelta();
            o.assertLessThan(norm(delta_matlab - deltay), tolerance);
            o.assertLessThan(norm(delta - deltay), tolerance); 
            %The test can fail, if the LDDi is badly conditioned.
        end
        
        function hesse_L = buildUpHesse(o, getLDD,n_var)
            % BUILDUPHESSE Builds up the big hesse matrix, should be used
            % for testing only, as this is quite inefficient. Note, that
            % this function works only, if every step has the same number
            % of additional constraints.
            hesse_L = zeros(n_var * (o.horizon +1) - o.n_contr - (n_var - 30) );
            for i = 1: o.horizon
                hesse_L( (i-1) * n_var +1:(i-1) * n_var +13  , (i-1) * n_var + 13+1: (i-1) * n_var + 26 ) = -eye(13);
                hesse_L( (i-1) * n_var + 13+1: (i-1) * n_var + 26 , (i-1) * n_var +1:(i-1) * n_var +13) = -eye(13);
                hesse_L( (i-1) * n_var +13+1 : (i-1)*n_var +13 + n_var, (i-1)*n_var +13 +1 : (i-1)*n_var +13 + n_var) = getLDD(i);
            end
            i = o.horizon +1 ;
            LDDi = getLDD(i);
            hesse_L( (i-1) * n_var +1:(i-1) * n_var +13  , (i-1) * n_var + 13+1: (i-1) * n_var + 26 ) = -eye(13);
            hesse_L( (i-1) * n_var + 13+1: (i-1) * n_var + 26 , (i-1) * n_var +1:(i-1) * n_var +13) = -eye(13);
            hesse_L( (i-1) * n_var +13+1 : end, (i-1)*n_var +13 +1 : end) = LDDi(1:13,1:13);
        end
        
        function delta = assembleDelta(o)
            delta = [];
            for i = 1: o.horizon
                y_t = [o.delta_lambda{i}; o.delta_s{i}; o.delta_q{i} ; o.delta_mu{i}];
                delta(end +1 : end + length(y_t)) = y_t;
            end
            i = o.horizon + 1;
            %Don't forget horizon +1
            y_t = [o.delta_lambda{i}; o.delta_s{i}];
            delta(end+1 : end + length(y_t)) = y_t;
            delta = delta';
        end
    end
end