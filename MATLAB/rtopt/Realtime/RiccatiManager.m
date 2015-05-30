classdef RiccatiManager < handle
    % RICCATISTEP Summary of this function goes here
    %   Detailed explanation goes here
    
    properties(Access=private)
        m ; %number of unknowns
        A ; %matrix storing parameters corresponding to smaller indices
        b ; %vector storing summands]
        
    end
    
    properties(Dependent)
        sol;
        solved;
    end
    
    methods
        function ricStep = RiccatiManager(m)
            obj.m = m;
            obj.A = zeros(m);
            obj.b = zeros(m,1);
            obj.sol = zeros(m,1);
            obj.solved = false;
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
            if(~obj.solved)
                obj.sol = obj.b;
                
                for i = 1:obj.m
                    obj.sol(i) = obj.A(i,:) * obj.sol;
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
        
        function sol = get.sol(obj)
            if obj.solved
                sol = obj.sol;
            else
                error('Recursion is not resolved yet');
            end
        end
    end
end