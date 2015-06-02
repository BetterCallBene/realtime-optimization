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
                end
            elseif(nargin == 1)
                obj.horizon = vargin{1};
                
                %Initialize storage
                obj.A = cell(horizon,1);
                obj.B = cell(horizon,1);
                obj.M = cell(horizon,1);
                obj.P = cell(horizon+1,1);
                obj.Q = cell(horizon+1,1);
                obj.R = cell(horizon,1);
                obj.nabla_s_star = cell(horizon+1, 1);
                obj.nabla_lambda = cell(horizon+1,1);
                obj.nabla_q = cell(horizon ,1);
                obj.delta = zeros((2*13+4)*horizon + 2*13, 1);
            else
                error('wrong number of inputs');
            end
        end
        
        function doStep(obj,i,LDD_i, LD_i)
            LD_i = -LD_i;
            obj.Q{i} = LDD_i(1:13,1:13);
            obj.nabla_lambda{i} = LD_i(1:13);
            obj.nabla_q{i} = LD_i(2*13+1:2*13+4);
            if(i == obj.horizon)
                obj.P{ i} = obj.Q{i};
                obj.nabla_s_star{i} = LD_i(13+1:2*13);
                
                
            else
                obj.M{i} = LDD_i(1:13,14:17);
                obj.R{i} = LDD_i(14:17, 14:17);
                
                obj.A{i} = LDD_i(18:end, 1:13);
                obj.B{i} = LDD_i(18:end, 14:17);
                
                obj.P{i} = obj.Q{i} + (obj.A{i}' * obj.P{i+1} * obj.A{i}) - ...
                    (obj.M{i} + obj.A{i}' * obj.P{i+1} * obj.B{i}) * ...
                    (obj.R{i} + obj.B{i}' * obj.P{i+1}  * obj.B{i}) \ ...
                    (obj.M{i}'  + obj.B{i}' * obj.P{i+1} * obj.A{i});
                
                obj.nabla_s_star{i} = LD_i(13+1:2*13) + obj.A{i}' * ...
                    obj.P{i+1} * obj.nabla_lambda{i+1} + obj.A{i}' *  ...
                    obj.nabla_s_star{i+1} -  ( obj.R{i} +  obj.B{i}' * ...
                    obj.P{i+1} *   obj.B{i}) ...
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
                obj.delta( (i-1) * 30 + 2*13 + 1 : i * 30 ) = (obj.R{i} + obj.B{i}' * obj.P{i+1} * obj.B{i}) \ ( obj.nabla_q{i} + obj.B{i}' * obj.P{i+1} * obj.nabla_lambda{i+1} + obj.B{i}' 
       
            
        end
        
    end
    
