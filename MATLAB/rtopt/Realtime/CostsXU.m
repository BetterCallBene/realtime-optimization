classdef CostsXU < Costs
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        alpha; % weight of the cost function for the control term
        beta;  % weight of the cost function for the state term
         cam_pos = [ 2; 0; 5];
    end
    
    methods
        % constructor
        function cC = CostsXU(varargin)
            cC@Costs(42,1337,pi); %Just some dummy variables
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
                % constructor based on two input values
                % a classDyn element and a classOCPparam element
            elseif (nargin >= 1)
                if (isa(varargin{1}, 'Dyn'))
                    cC.dyn = varargin{1};
                else
                    error('wrong class type for problem parameter');
                end
                
                if (nargin == 3)
                    cC.alpha = varargin{2};
                    cC.beta = varargin{3};
                    
                else
                    error('wrong number of inputs');
                end
            end
        end
        
        function c_val = get_cost(obj)
            % compute the cost value with penalty term
            control = obj.dyn.contr;
            state   = obj.dyn.state;
           
            one_vec = ones(1,size(state,2));
            
            alpha_   = obj.alpha;
            beta_    = obj.beta;
            
            mesh    = obj.dyn.environment.mesh;
            c_val   = alpha_* 0.5*sum((control(:,1:end-1).^2)*mesh') + ...
                beta_ * 0.5*sum(((state(1:3,1:end) - kron(one_vec,obj.cam_pos)).^2)*one_vec');
        end
        
        function cD_val = get_costD(obj)
            % compute the Jacobian of the cost value with penalty term
            [n_tp, n_state, n_contr] = getParams(obj);
            
            control     = obj.dyn.contr;
            mesh        = obj.dyn.environment.mesh;
            state   = obj.dyn.state;
            
            alpha_   = obj.alpha;
            beta_    = obj.beta;
            
            rind        = zeros(1,(n_tp-1)*n_contr);
            vind        = zeros(1,(n_tp-1)*n_contr);
            
            rindx       = zeros(1,(n_tp)*3);
            vindx       = zeros(1,(n_tp)*3);
            for i = 1:n_tp-1
                rind((i-1)*n_contr+1:i*n_contr) = ...
                    (i-1)*(n_state+n_contr)+n_state + 1:...
                    i*(n_state+n_contr);
                vind((i-1)*n_contr+1:i*n_contr) = ...
                    alpha_ * mesh(i)*control(:,i)';
                
                rindx((i-1)*3+1:i*3) = (i-1)*(n_state+n_contr) + 1:...
                    (i-1)*(n_state+n_contr)+3;
                vindx((i-1)*3+1:i*3) = beta_ * (state(1:3,i) - obj.cam_pos);
            end
            
            i = n_tp;
            rindx((i-1)*3+1:i*3) = (i-1)*(n_state+n_contr) + 1:...
                (i-1)*(n_state+n_contr)+3;
            vindx((i-1)*3+1:i*3) = beta_ * (state(1:3,i) - obj.cam_pos);
            
            cD_val        = sparse(rind,ones(1,(n_tp-1)*n_contr),vind,...
                n_tp*(n_state+n_contr),1) + ...
                sparse(rindx,ones(1,n_tp*3),vindx,n_tp*(n_state+n_contr),1);
        end
        
        function cDD_val = get_costDD(obj)
            %compute the Hessian of the cost value with penalty term
            [n_tp, n_state, n_contr] = getParams(obj);
            
            mesh        = obj.dyn.environment.mesh;
            
            alpha_   = obj.alpha;
            beta_    = obj.beta;
            
            rind        = zeros(1,(n_tp-1)*n_contr);
            vind        = zeros(1,(n_tp-1)*n_contr);
            
            rindx       = zeros(1,(n_tp)*3);
            vindx       = zeros(1,(n_tp)*3);
            for i = 1:n_tp-1
                rind((i-1)*n_contr+1:i*n_contr) = ...
                    (i-1)*(n_state+n_contr)+n_state + 1:...
                    i*(n_state+n_contr);
                vind((i-1)*n_contr+1:i*n_contr) = ...
                    alpha_ * mesh(i)*ones(1,n_contr);
                
                rindx((i-1)*3+1:i*3) = (i-1)*(n_state+n_contr) + 1:...
                    (i-1)*(n_state+n_contr)+3;
                vindx((i-1)*3+1:i*3) = beta_ * ones(1,3);
            end
            
            i = n_tp;
            rindx((i-1)*3+1:i*3) = (i-1)*(n_state+n_contr) + 1:...
                (i-1)*(n_state+n_contr)+3;
            vindx((i-1)*3+1:i*3) = beta_ * ones(1,3);
            
            cDD_val        = {sparse(rind,rind,vind,...
                n_tp*(n_state+n_contr),n_tp*(n_state+n_contr)) + ...
                sparse(rindx,rindx,vindx,n_tp*(n_state+n_contr),n_tp*(n_state+n_contr))};
        end
        
        function cDD_val = get_costDD_approx(obj,t)
            %compute the approximated Hessian of the cost value
            [n_tp, n_state, n_contr] = getParams(obj);
            
            mesh        = obj.dyn.environment.mesh;
            
            alpha_ = obj.alpha;
            beta_    = obj.beta;
            if t >= n_tp
                CostD_t = [beta_*ones(3,1);sparse(n_state-3,1); zeros(n_contr,1)];
            else
                CostD_t = [beta_*ones(3,1);sparse(n_state-3,1); alpha_ * mesh(t)*ones(n_contr,1)];
            end
            cDD_val = CostD_t * CostD_t' ;
            
        end
        
    end
    
    methods
        
        function setupTest(o,horizon)
            env = Environment();
            env.wind = @(s_t, t)  s_t + [rand(3,1); zeros(10,1)];
%             env.setUniformMesh(uint16(horizon));
            env.setUniformMesh1(horizon +1,1); 
            cQ = Quadrocopter();
            
            tol = 1e-3;
            opts = odeset('RelTol',tol,'AbsTol',0.1*tol);
            cFE = ode15sM(opts);
            
            % Initialisierung der Dynamik
            cBQD = BasisQDyn(cQ, env,cFE);
            cBQD.vec = rand(cQ.n_var * (horizon +1),1);
            
            % CostsXU parametrisieren     
            o.dyn = cBQD;
            o.vec = o.dyn.vec;
            o.alpha = 2;
            o.beta = 3;
            
        end
    end
    
    
    
end

