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
        horizon ;
        delta;
        
        n_lambda;
        n_state;
        n_contr;
        n_var; 
        
    end
    
    methods
        function ricM = RiccatiManager_woConstr(varargin)
            % RICCATIMANAGER_WOCONSTR Takes as inpt the horizon length and
            % the used model. If not specified, the defauls are 20 and for
            % the model a instance of Quadrocopter is used.
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                else
                    obj.horizon = 20;
                    obj.robot = Quadrocopter();
                end
            elseif(nargin == 1)
                obj.horizon = varargin{1};
                obj.robot = Quadrocopter();
            elseif(nargin == 2)
                obj.horizon = varargin{1};
                if (isa(varargin{2},'Model'))
                obj.robot = varargin{2};
                end
            else
                error('wrong number of inputs');
            end
           
            %Initialize storage
            obj.n_lambda = obj.robot.n_state;
            obj.n_state = obj.robot.n_state;
            obj.n_contr = obj.robot.n_contr;
            obj.n_var = obj.n_lambda + obj.n_state + obj.n_contr;
            
            obj.A = cell(obj.horizon,1);
            obj.B = cell(obj.horizon,1);
            obj.M = cell(obj.horizon,1);
            obj.P = cell(obj.horizon+1,1);
            obj.Q = cell(obj.horizon+1,1);
            obj.R = cell(obj.horizon,1);
            obj.nabla_s_star = cell(obj.horizon+1, 1);
            obj.nabla_lambda = cell(obj.horizon+1,1);
            obj.nabla_q = cell(obj.horizon ,1);
            obj.delta = zeros( obj.n_var * (obj.horizon +1 ) - obj.n_contr, 1);
            
            
            
        end
        
        function doStep(o,i,LDD_i, LD_i)
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
            
            %solve [delta lambda ; delta s]
            if ( i ==1 )
                o.delta(1: o.n_lambda + o.n_state) = ...
                    [-o.P{i} , -eye(o.n_lambda) ; -eye(o.n_state) , zeros(o.n_state)] * ...
                    [ o.nabla_lambda{i} ; o.nabla_s_star{i} ];
            else
                o.delta((i-1) * 30 + 1 : (i-1) * 30 + 2*13) = [-o.P{i} , -eye(o.n_lambda) ; -eye(o.n_state) , zeros(o.n_state)] * [o.nabla_lambda{i} - o.A{i-1} * o.delta( (i-2) * 30 + 13 +1 : (i-2) * 30 +2*13) + o.B{i-1} * o.delta( (i-2) *30 + 2*13 +1: (i-1)*30) ; o.nabla_s_star{i} ];
            end
          
            %solve for delta q
            if( i ~= o.horizon +1 )
            o.delta( (i-1) * 30 + 2*13 + 1 : i * 30 ) = ...
                (o.R{i} + o.B{i}' * o.P{i+1} * o.B{i}) \ ...
                ( ...
                    o.nabla_q{i} + ...
                    o.B{i}' * o.P{i+1} * o.nabla_lambda{i+1} + ...
                    o.B{i}' * o.nabla_s_star{i+1} - ...
                    ( o.M{i}' + o.B{i}' * o.P{i+1} * o.A{i} ) * o.delta( (i-1) * 30 + 13 + 1: (i-1)*30 + 2*13) ...
                );
            
            % nach dem ersten Schritt kann dann u = q_t + deltaq_t
            % �bergeben werden um zu steuern
            end
        end
    end
    
    methods(Test)
        
%         function testRiccati_0(obj)
% %             for i = 1:100
%             obj.doTestWithHorizon(0);
% %             end
%         end
        
        function testRiccati_1(obj)
            obj.doTestWithHorizon(1);
        end
        
%         function testRiccati_100(obj)
%             obj.doTestWithHorizon(100);
%         end
    end
    
    methods
        function initialize(obj)
            %Initialize storage
            obj.n_lambda = obj.robot.n_state;
            obj.n_state = obj.robot.n_state;
            obj.n_contr = obj.robot.n_contr;
            obj.n_var = obj.n_lambda + obj.n_state + obj.n_contr;
            
            obj.A = cell(obj.horizon,1);
            obj.B = cell(obj.horizon,1);
            obj.M = cell(obj.horizon,1);
            obj.P = cell(obj.horizon+1,1);
            obj.Q = cell(obj.horizon+1,1);
            obj.R = cell(obj.horizon,1);
            obj.nabla_s_star = cell(obj.horizon+1, 1);
            obj.nabla_lambda = cell(obj.horizon+1,1);
            obj.nabla_q = cell(obj.horizon ,1);
            obj.delta = zeros( obj.n_var * (obj.horizon +1 ) - obj.n_contr, 1);
        
        end
        
        function doTestWithHorizon(obj, hori)
            obj.horizon = hori;
            obj.robot  = Quadrocopter();
            obj.initialize();
            
            LDDi = zeros(30);
            %Find a invertible matrix
            while (det(LDDi) == 0  || abs(det(LDDi)) > 1e2)
                LDDi = 10 * eye(30) +  rand(30);
                LDDi(18:end, 18:end) = zeros(13);
                %LDDi must be symmetric
                LDDi = .5* (LDDi + LDDi');
            end
            
            getLDD = @(i) LDDi;
            hesse_L = obj.buildUpHesse(getLDD);
            
            deltay = 10 * rand((obj.horizon+1) * 30 - 4,1);
            grad_L = hesse_L * deltay;
            
            %Perform reccati recursion
            for i = obj.horizon+1:-1:1
                if i == obj.horizon +1
                    obj.doStep(i,LDDi, grad_L( (i-1) * 30 +1 : end ,1  ) );
                else
                    obj.doStep(i,LDDi, grad_L( (i-1) * 30 +1 : i * 30,1  ));
                end
            end
            
            for i = 1:obj.horizon +1
                obj.solveStep(i);
            end
            
            obj.assertLessThan(norm(obj.delta - deltay), obj.tol);
        end
        
        function hesse_L = buildUpHesse(obj, getLDD)
            hesse_L = zeros(30 * obj.horizon +26);
            for i = 1: obj.horizon
                hesse_L( (i-1) * 30 +1:(i-1) * 30 +13  , (i-1) * 30 + 13+1: (i-1) * 30 + 26 ) = -eye(13);
                hesse_L( (i-1) * 30 + 13+1: (i-1) * 30 + 26 , (i-1) * 30 +1:(i-1) * 30 +13) = -eye(13);
                hesse_L( (i-1) * 30 +13+1 : (i-1)*30 +13 + 30, (i-1)*30 +13 +1 : (i-1)*30 +13 + 30) = getLDD(i);
            end
            
            i = obj.horizon +1 ;
            LDDi = getLDD(i);
            hesse_L( (i-1) * 30 +1:(i-1) * 30 +13  , (i-1) * 30 + 13+1: (i-1) * 30 + 26 ) = -eye(13);
            hesse_L( (i-1) * 30 + 13+1: (i-1) * 30 + 26 , (i-1) * 30 +1:(i-1) * 30 +13) = -eye(13);
            hesse_L( (i-1) * 30 +13+1 : end, (i-1)*30 +13 +1 : end) = LDDi(1:13,1:13);
        end
    end
end
