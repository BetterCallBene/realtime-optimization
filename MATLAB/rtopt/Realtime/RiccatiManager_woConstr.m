classdef RiccatiManager_woConstr < TestEnv
    % RICCATIMANAGER_WOCONSTR Implements the Riccati algorithm without the constraints
    
    
    properties(Access=private)
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
    end
    
    methods
        function ricM = RiccatiManager_woConstr(varargin)
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                else
                    obj.horizon = 20;
                end
            elseif(nargin == 1)
                obj.horizon = varargin{1};
            else
                error('wrong number of inputs');
            end
            %Initialize storage
            obj.A = cell(obj.horizon,1);
            obj.B = cell(obj.horizon,1);
            obj.M = cell(obj.horizon,1);
            obj.P = cell(obj.horizon+1,1);
            obj.Q = cell(obj.horizon+1,1);
            obj.R = cell(obj.horizon,1);
            obj.nabla_s_star = cell(obj.horizon+1, 1);
            obj.nabla_lambda = cell(obj.horizon+1,1);
            obj.nabla_q = cell(obj.horizon ,1);
            obj.delta = zeros( (2*13+4)*obj.horizon + 2*13, 1);
        end
        
        function doStep(obj,i,LDD_i, LD_i)
            LD_i = -LD_i;
            obj.Q{i} = LDD_i(1:13,1:13);
            obj.nabla_lambda{i} = LD_i(1:13);
            
            
            if(i == obj.horizon+1)
                obj.P{i} = obj.Q{i};
                obj.nabla_s_star{i} = LD_i(13+1:2*13);
                
            else
                obj.nabla_q{i} = LD_i(2*13+1:end);
%             obj.nabla_q{i} = LD_i(2*13+1:2*13+4);
                obj.M{i} = LDD_i(1:13,14:17);
                obj.R{i} = LDD_i(14:17, 14:17);
                
                obj.A{i} = LDD_i(18:end, 1:13);
                obj.B{i} = LDD_i(18:end, 14:17);
                
                obj.P{i} = obj.Q{i} + (obj.A{i}' * obj.P{i+1} * obj.A{i}) - ...
                    (obj.M{i} + obj.A{i}' * obj.P{i+1} * obj.B{i}) * ...
                    ((obj.R{i} + obj.B{i}' * obj.P{i+1}  * obj.B{i}) \ ...
                    (obj.M{i}'  + obj.B{i}' * obj.P{i+1} * obj.A{i}));
                
                obj.nabla_s_star{i} = LD_i(13+1:2*13) + ...
                    obj.A{i}' * obj.P{i+1} * obj.nabla_lambda{i+1} + ...
                    obj.A{i}' * obj.nabla_s_star{i+1} - ...
                    ( obj.R{i} +  obj.B{i}' * obj.P{i+1} *   obj.B{i}) ...
                    \ ( obj.nabla_q{i} +  obj.B{i}' *   obj.P{i+1} *  obj.nabla_lambda{i+1} +  obj.B{i}' *  obj.nabla_s_star{i+1});
                
            end
        end
        
        function solveStep(obj, i)
            
            %solve delta lambda
            if ( i ==1 )
                obj.delta(1:2*13) = [-obj.P{i} , -spones(13) ; -spones(13) , spzeros(13)] * [ obj.nabla_lambda{i} ; obj.nabla_s_star{i} ];
            else
                obj.delta((i-1) * 30 + 1 : (i-1) * 30 + 2*13) = [-obj.P{i} , -spones(13) ; -spones(13) , spzeros(13)] * [obj.nabla_lambda{i} - obj.A{i-1} * obj.delta( (i-2) * 30 + 13 +1 : (i-2) * 30 +2*13) + obj.B{i-1} * obj.delta( (i-2) *30 + 2*13 +1: (i-1)*30) ; obj.nabla_s_star{i} ];
            end
            obj.delta( (i-1) * 30 + 2*13 + 1 : i * 30 ) = (obj.R{i} + obj.B{i}' * obj.P{i+1} * obj.B{i}) \ ( obj.nabla_q{i} + obj.B{i}' * obj.P{i+1} * obj.nabla_lambda{i+1} + obj.B{i}' * obj.nabla_s_star{i+1} - ( obj.M{i}' + obj.B{i}' * obj.P{i+1} * obj.A{i}) * obj.delta( (i-1) * 30 + 13 + 1: (i-1)*30 + 2*13));
            
            % nach dem ersten Schritt kann dann u = q_t + deltaq_t
            % �bergeben werden um zu steuern
            
        end
        
    end
    
    methods(Test)
        
        function testRiccati_ow(obj)
            obj.horizon = 20;
            %Initialize storage
            obj.A = cell(obj.horizon,1);
            obj.B = cell(obj.horizon,1);
            obj.M = cell(obj.horizon,1);
            obj.P = cell(obj.horizon+1,1);
            obj.Q = cell(obj.horizon+1,1);
            obj.R = cell(obj.horizon,1);
            obj.nabla_s_star = cell(obj.horizon+1, 1);
            obj.nabla_lambda = cell(obj.horizon+1,1);
            obj.nabla_q = cell(obj.horizon ,1);
            obj.delta = zeros( (2*13+4)*obj.horizon + 2*13, 1);
            
            
            LDDi = zeros(30);
            %Find a invertible matrix
            while (det(LDDi) == 0  || abs(det(LDDi)) > 1e2)
                LDDi = 10 * eye(30) +  rand(30);
                LDDi(18:end, 18:end) = zeros(13);
            end
            
            hesse_L = zeros(30 * obj.horizon +26);
            for i = 1: obj.horizon
                hesse_L( (i-1) * 30 +1:(i-1) * 30 +13  , (i-1) * 30 + 13+1: (i-1) * 30 + 26 ) = -eye(13);
                hesse_L( (i-1) * 30 + 13+1: (i-1) * 30 + 26 , (i-1) * 30 +1:(i-1) * 30 +13) = -eye(13);
                hesse_L( (i-1) * 30 +13+1 : (i-1)*30 +13 + 30, (i-1)*30 +13 +1 : (i-1)*30 +13 + 30) = LDDi;
            end
            
            i = obj.horizon +1 ;
            hesse_L( (i-1) * 30 +1:(i-1) * 30 +13  , (i-1) * 30 + 13+1: (i-1) * 30 + 26 ) = -eye(13);
            hesse_L( (i-1) * 30 + 13+1: (i-1) * 30 + 26 , (i-1) * 30 +1:(i-1) * 30 +13) = -eye(13);
            hesse_L( (i-1) * 30 +13+1 : end, (i-1)*30 +13 +1 : end) = LDDi(1:13,1:13);
            
            
            deltay = 10 * rand(626,1);
            grad_L = - hesse_L * deltay;
            
            
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
        
    end
end
