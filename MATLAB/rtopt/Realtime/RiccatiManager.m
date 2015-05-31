classdef RiccatiManager <  TestEnv
    % RICCATIMANAGER Central class to perform a riccati solution
    
    properties(Access=private)
        m ; %number of unknowns
        A ; %matrix storing parameters corresponding to smaller indices
        b ; %vector storing summands]
        solved;
        sol;
    end
    
    properties(Dependent)
        
        
    end
    
    methods
        function ricStep = RiccatiManager(varargin)
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
            elseif(nargin == 1)
                obj.m = varargin{1};
                obj.A = zeros(varargin{1});
                obj.b = zeros(varargin{1},1);
                obj.sol = zeros(varargin{1},1);
                obj.solved = false;
            else
                error('wrong number of inputs');
            end
        end
        
        function doStep(obj,i, LDD_i, LD_i, lim_right)
            % DOSTEP This function expresses the i_th component of the unknown as a
            % sum of unknowns with smaller indices. Please ensure, that in LDD i is
            % the largest non zero element, not mentioned in the recursion yet.
            
            obj.A(i,:) = -LDD_i/LDD_i(i);
            obj.b(i) = LD_i / LDD_i(i);
            
            obj.A(i,i) = 0;
            
            % apply recursion
            for j = min(i+lim_right, obj.m):-1: i
                if obj.A(i,j) ~= 0
                    
                    obj.A(i,:) = obj.A(i,:) - (obj.A(i,j)* obj.A(j,:));
                    obj.b(i) = obj.b(i) - (obj.A(i,j) * obj.b(j));
                    
                    obj.A(i,j) = 0;
                end
            end
        end
        
        function res = resolveRecursion(obj)
            obj.A = obj.A + eye(obj.m);
            if(~obj.solved)
                obj.sol = obj.b;
                
                for i = 1:obj.m
                    obj.sol(i)
                    obj.sol(i) = obj.A(i,:) * obj.sol;
                    obj.sol(i)
                end
            end
            res = obj.sol;
            obj.solved = true;
        end
        
        function res = solveFirst_k_elements(obj, k)
            obj.sol = obj.b;
            for i = 1:min(k,obj.m)
                obj.sol(i) = obj.A(i,:) * obj.sol;
            end
            res = obj.sol(1:k);
        end
        
        function res = solveFrom_kp1_elements(obj, k)
            for i = k+1:obj.m
                obj.sol(i) = obj.A(i,:) * obj.sol;
            end
            res = obj.sol(k+1:end);
            obj.solved = true;
        end
        
        %         function sol = get.sol(obj)
        %             if obj.solved
        %                 sol = obj.sol;
        %             else
        %                 sol = obj.sol;
        %                 disp('Recursion is not resolved yet');
        %             end
        %         end
    end
    
    methods(Test)
        
        function testRiccati(obj)
            
            %Test algorithm with some random matrix A and vector b
            
            testm = 5;
            while(true)
                testA = rand(testm)*100+ ones(testm);
                if(det(testA) > 0.1)
                    break;
                end
            end
            obj.m = testm;
            obj.A = testA;
            testb = rand(testm,1)*100+ones(testm,1);
            obj.b = testb;
            obj.solved = false;
            
            
            %Calculate numerical solution
            numX = testA\testb;
            
            %Calculate riccati solution
            for i = testm:-1:1
                obj.doStep(i, testA(i,:), testb(i), obj.m);
            end
            
            ricX = obj.resolveRecursion();
            
            obj.assertLessThan(norm(ricX - numX), obj.tol);
        end
        
    end
end