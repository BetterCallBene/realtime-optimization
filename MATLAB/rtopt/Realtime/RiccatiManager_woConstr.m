classdef RiccatiManager_woConstr < TestEnv
    % RICCATIMANAGER_WOCONSTR Implements the Riccati algorithm without
    % additional constraints.
    
    properties(Access=private)
        robot;
        A;
        B;
        M;
        P;
        Q;
        R;
        nabla_s_star;
        nabla_lambda;
        nabla_q;
        horizon;
                
               
        n_lambda;
        n_state;
        n_contr;
        n_var;
    end
    
    properties
        delta;
        
        delta_lambda;
        delta_s;
        delta_q;
    end
    
    methods
        function o = RiccatiManager_woConstr(varargin)
            % RICCATIMANAGER_WOCONSTR Takes as inpt the horizon length and
            % the used model. If not specified, the defauls are 20 and for
            % the model a instance of Quadrocopter.
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
            
            %Initialize storage
            o.initialize();
        end
        
        function doStep(o,i,LDD_i, LD_i)
            % DOSTEP Performs one step of Riccati Recursion Part 1 
            o.Q{i} = LDD_i(1:o.n_state,1:o.n_state);
            o.nabla_lambda{i} = LD_i(1:o.n_lambda);
            
            if(i == o.horizon+1)
                o.P{i} = o.Q{i};
                o.nabla_s_star{i} = LD_i(o.n_lambda+1:o.n_lambda + o.n_state);
                
            else
                o.nabla_q{i} = LD_i(o.n_lambda + o.n_state +1  : o.n_var);
                o.M{i} = LDD_i(1:o.n_state,o.n_state +1 :o.n_state + o.n_contr);
                o.R{i} = LDD_i(o.n_state +1 :o.n_state + o.n_contr, o.n_state+1 :o.n_state + o.n_contr);
                
                o.A{i} = LDD_i(o.n_state + o.n_contr +1: o.n_var, 1:o.n_state);
                o.B{i} = LDD_i(o.n_state + o.n_contr +1: o.n_var, o.n_state +1 : o.n_state + o.n_contr);
                
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
                
            end
        end
        
        function solveStep(o, i)
            % SOLVESTEP Resolves one step of Riccati Recursion Part 2
            
            %solve [delta lambda ; delta s]
            if ( i ==1 )
                
                o.delta_lambda{i} = -o.P{i}* o.nabla_lambda{i} - o.nabla_s_star{i};
                o.delta_s{i} = -o.nabla_lambda{i};
                
                o.delta(1: o.n_lambda + o.n_state) = [o.delta_lambda{i}; o.delta_s{i}];
                
            else
                var = (o.A{i-1} * o.delta_s{i-1} + o.B{i-1} * o.delta_q{i-1});
                
                % solve delta lambda
                o.delta_lambda{i} = -o.P{i}* o.nabla_lambda{i} - o.nabla_s_star{i} + o.P{i} * var ;
                % solve delta s
                o.delta_s{i} = - o.nabla_lambda{i}  + var;
                
                o.delta((i-1) * 30 + 1 : (i-1) * 30 + 2*13) = [o.delta_lambda{i} ; o.delta_s{i}];
                
            end
            
            %solve for delta q
            if( i ~= o.horizon +1 )
                o.delta_q{i} = (o.R{i} + o.B{i}' * o.P{i+1} * o.B{i}) \ ...
                    ( ...
                    o.nabla_q{i} + ...
                    o.B{i}' * o.P{i+1} * o.nabla_lambda{i+1} + ...
                    o.B{i}' * o.nabla_s_star{i+1} - ...
                    ( o.M{i}' + o.B{i}' * o.P{i+1} * o.A{i} ) * o.delta_s{i} );
                
                
                o.delta( (i-1) * 30 + 2*13 + 1 : i * 30 ) = o.delta_q{i};
            end
        end
    end
    
    methods(Test)
        
        function testRiccati_0(o)
            o.doTestWithHorizon(0);
        end
        
        function testRiccati_1(o)
            o.doTestWithHorizon(1);
        end
        
        function testRiccati_20(o)
            o.doTestWithHorizon(20);
        end
        
        function testRiccati_100(o)
            o.doTestWithHorizon(100);
        end
    end
    
    methods
        function initialize(o)
            %Initialize storage
            o.n_lambda = o.robot.n_state;
            o.n_state = o.robot.n_state;
            o.n_contr = o.robot.n_contr;
            o.n_var = o.n_lambda + o.n_state + o.n_contr;
            
            o.A = cell(o.horizon,1);
            o.B = cell(o.horizon,1);
            o.M = cell(o.horizon,1);
            o.P = cell(o.horizon+1,1);
            o.Q = cell(o.horizon+1,1);
            o.R = cell(o.horizon,1);
            o.nabla_s_star = cell(o.horizon+1, 1);
            o.nabla_lambda = cell(o.horizon+1,1);
            o.nabla_q = cell(o.horizon ,1);
            o.delta = zeros( o.n_var * (o.horizon +1 ) - o.n_contr, 1);
            
        end
        
        function doTestWithHorizon(o, hori)
            o.horizon = hori;
            o.robot  = Quadrocopter();
            o.initialize();
            
            LDDi = zeros(30);
            %Find a invertible matrix
            while (det(LDDi) == 0  || abs(det(LDDi)) > 1e2)
                LDDi = 10 * eye(30) +  rand(30);
                LDDi(18:end, 18:end) = zeros(13);
                %LDDi must be symmetric
                LDDi = .5* (LDDi + LDDi');
            end
            
            getLDD = @(i) LDDi;
            hesse_L = o.buildUpHesse(getLDD);
            
            deltay = 10 * rand((o.horizon+1) * 30 - 4,1);
            grad_L = hesse_L * deltay;
            
            disp(' ');
            
            tic;
            %Perform reccati recursion
            for i = o.horizon+1:-1:1
                if i == o.horizon +1
                    o.doStep(i,LDDi, grad_L( (i-1) * 30 +1 : end ,1  ) );
                else
                    o.doStep(i,LDDi, grad_L( (i-1) * 30 +1 : i * 30,1  ));
                end
            end
            
            for i = 1:o.horizon +1
                o.solveStep(i);
            end
            timeRiccati = toc;
            disp(strcat('Time for Riccati with Horizon',{' '}, int2str(hori), ':'));
            disp(num2str(timeRiccati));
            
            tic;
            hesse_L = o.buildUpHesse(getLDD);
            delta_matlab = hesse_L \ grad_L;
            timeMatlab = toc;
            disp(strcat('Time for \ with Horizon',{' '}, int2str(hori), ': '));
            disp(num2str(timeMatlab));
            
            o.assertLessThan(norm(delta_matlab - deltay), o.tol * 10);
            o.assertLessThan(norm(o.delta - deltay), o.tol * 10);
        end
        
        function hesse_L = buildUpHesse(o, getLDD)
            % BUILDUPHESSE Builds up the big hesse matrix, should be used
            % for testing only, as this is quite inefficient
            hesse_L = zeros(o.n_var * (o.horizon +1) - o.n_contr);
            for i = 1: o.horizon
                hesse_L( (i-1) * 30 +1:(i-1) * 30 +13  , (i-1) * 30 + 13+1: (i-1) * 30 + 26 ) = -eye(13);
                hesse_L( (i-1) * 30 + 13+1: (i-1) * 30 + 26 , (i-1) * 30 +1:(i-1) * 30 +13) = -eye(13);
                hesse_L( (i-1) * 30 +13+1 : (i-1)*30 +13 + 30, (i-1)*30 +13 +1 : (i-1)*30 +13 + 30) = getLDD(i);
            end
            
            i = o.horizon +1 ;
            LDDi = getLDD(i);
            hesse_L( (i-1) * 30 +1:(i-1) * 30 +13  , (i-1) * 30 + 13+1: (i-1) * 30 + 26 ) = -eye(13);
            hesse_L( (i-1) * 30 + 13+1: (i-1) * 30 + 26 , (i-1) * 30 +1:(i-1) * 30 +13) = -eye(13);
            hesse_L( (i-1) * 30 +13+1 : end, (i-1)*30 +13 +1 : end) = LDDi(1:13,1:13);
        end
        
        function hesse_L = buildUpHesse_activ(o, getLDD)
            % BUILDUPHESSE Builds up the big hesse matrix with activ constraints, should be used
            % for testing only, as this is quite inefficient
            
            LDD1 = getLDD(1);
            [m , n] = size(LDD1);
            hesse_L = [sparse(13,13), -speye(13), sparse(13,n-13);...
                       [ -speye(13); sparse(m-13,13)], LDD1];
            
            for i = 2: o.horizon
                LDDi = getLDD(i);
                [m , n] = size(hesse_L);
                [k , l] = size(LDDi);
                
                SPE = [ sparse(13,n-13),-speye(13);...
                        sparse(k-13,n)];
                hesse_L = [hesse_L, SPE';
                           SPE, LDDi];
            end
            
            i = o.horizon +1 ;
            LDDi = getLDD(i);
            [m, n] = size(hesse_L);
            SPE = [ sparse(13,n-13),-speye(13)];
            hesse_L = [hesse_L, SPE';
                       SPE, LDDi(1:13,1:13)];
        end
        
        function grad_L = buildUpGradient(o, getLD)
            % BUILDUPGRADIENT Builds the huge gradient, should be used for
            % testing only.
            grad_L = zeros(o.n_var * (o.horizon +1) - o.n_contr,1);
            for i = 1:o.horizon
                grad_L( (i-1) * o.n_var +1 : i * o.n_var) = getLD(i);
            end
            i = o.horizon +1;
            LDi = getLD(i);
            grad_L( (i-1 ) * o.horizon +1 : end) = LDi(1: o.n_state+o.n_lambda);
        end
        
        function grad_L = buildUpGradient_activ(o, getLD)
            % BUILDUPGRADIENT Builds the huge gradient with activ constraints, should be used for
            % testing only.
            grad_L = getLD(1);
            for i = 2:o.horizon
                grad_L = [grad_L , getLD(i)];
            end
            i = o.horizon +1;
            LDi = getLD(i);
            grad_L = [grad_L ,LDi(1: o.n_state+o.n_lambda)];
        end
    end
end
