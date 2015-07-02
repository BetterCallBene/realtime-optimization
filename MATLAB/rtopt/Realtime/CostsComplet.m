classdef CostsComplet < Costs
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        alpha; % weight of the cost function for the control term
        beta;  % weight of the cost function for the state term
        gamma; % weight of the cost function for the quadrions
        kappa; % weight of the cost function for the velocitys
        
        cam_pos 
        
        timepoint;
    end
    
    methods
        % constructor
        function cC = CostsComplet(varargin)
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
                    cC.cam_pos = @(timepoint) [ 2; 0; 5];
                else
                    error('wrong class type for problem parameter');
                end
                
                if (nargin == 5)
                    cC.alpha = varargin{2};
                    cC.beta = varargin{3};
                    cC.gamma = varargin{4};
                    cC.kappa = varargin{5};
                    
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
            gamma_   = obj.gamma;
            kappa_   = obj.kappa;
            
            mesh    = obj.dyn.environment.mesh;
            c_val   = alpha_* 0.5*sum((control(:,1:end-1).^2)*mesh') + ...
                beta_ * 0.5*sum(((state(1:3,1:end) - kron(one_vec,obj.cam_pos(obj.timepoint) )).^2)*one_vec') + ...
                gamma_ * 0.5*sum((state(4:7,1:end).^2)*one_vec') + ...
                kappa_ * 0.5*sum((state(8:end,1:end).^2)*one_vec');
        end
        
        function cD_val = get_costD(obj)
            % compute the Jacobian of the cost value with penalty term
            [n_tp, n_state, n_contr] = getParams(obj);
            
            control     = obj.dyn.contr;
            mesh        = obj.dyn.environment.mesh;
            state   = obj.dyn.state;
            
            
            alpha_   = obj.alpha;
            beta_    = obj.beta;
            gamma_   = obj.gamma;
            kappa_   = obj.kappa;
            
            rind        = zeros(1,(n_tp-1)*n_contr);
            vind        = zeros(1,(n_tp-1)*n_contr);
            
            rindx       = zeros(1,(n_tp)*n_state);
            vindx       = zeros(1,(n_tp)*n_state);
            for i = 1:n_tp-1
                rind((i-1)*n_contr+1:i*n_contr) = ...
                    (i-1)*(n_state+n_contr)+n_state + 1:...
                    i*(n_state+n_contr);
                vind((i-1)*n_contr+1:i*n_contr) = ...
                    alpha_ * mesh(i)*control(:,i)';
                
                rindx((i-1)*n_state+1:i*n_state) = (i-1)*(n_state+n_contr) + 1:...
                    (i-1)*(n_state+n_contr)+n_state;
                vindx((i-1)*n_state+1:(i-1)*n_state+3) = beta_ * (state(1:3,i) - obj.cam_pos(obj.timepoint));
                vindx((i-1)*n_state+4:(i-1)*n_state+7) = gamma_* state(4:7,i);
                vindx((i-1)*n_state+8:i*n_state)       = kappa_* state(8:end,i);
            end
            
            i = n_tp;
            rindx((i-1)*n_state+1:i*n_state) = (i-1)*(n_state+n_contr) + 1:...
                    (i-1)*(n_state+n_contr)+n_state;
            vindx((i-1)*n_state+1:(i-1)*n_state+3) = beta_ * (state(1:3,i) - obj.cam_pos(obj.timepoint));
            vindx((i-1)*n_state+4:(i-1)*n_state+7) = gamma_* state(4:7,i);
            vindx((i-1)*n_state+8:i*n_state)       = kappa_* state(8:end,i);
            
            cD_val = sparse(rind,ones(1,(n_tp-1)*n_contr),vind,...
                n_tp*(n_state+n_contr),1) + ...
                sparse(rindx,ones(1,n_tp*n_state),vindx,n_tp*(n_state+n_contr),1);
        end
        
        function cDD_val = get_costDD(obj)
            %compute the Hessian of the cost value with penalty term
            [n_tp, n_state, n_contr] = getParams(obj);
            
            mesh        = obj.dyn.environment.mesh;
            
            alpha_   = obj.alpha;
            beta_    = obj.beta;
            gamma_   = obj.gamma;
            kappa_   = obj.kappa;
            
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
                
                rindx((i-1)*n_state+1:i*n_state) = (i-1)*(n_state+n_contr) + 1:...
                    (i-1)*(n_state+n_contr)+n_state;
                vindx((i-1)*n_state+1:(i-1)*n_state+3) = beta_ * ones(1,3);
                vindx((i-1)*n_state+4:(i-1)*n_state+7) = gamma_* ones(1,4);
                vindx((i-1)*n_state+8:i*n_state)       = kappa_* ones(1,6);
            end
            
            i = n_tp;
            rindx((i-1)*n_state+1:i*n_state) = (i-1)*(n_state+n_contr) + 1:...
                    (i-1)*(n_state+n_contr)+n_state;
            vindx((i-1)*n_state+1:(i-1)*n_state+3) = beta_ * ones(1,3);
            vindx((i-1)*n_state+4:(i-1)*n_state+7) = gamma_* ones(1,4);
            vindx((i-1)*n_state+8:i*n_state)       = kappa_* ones(1,6);
            
            cDD_val        = {sparse(rind,rind,vind,...
                n_tp*(n_state+n_contr),n_tp*(n_state+n_contr)) + ...
                sparse(rindx,rindx,vindx,n_tp*(n_state+n_contr),n_tp*(n_state+n_contr))};
        end
        
        function cDD_val = get_costDD_approx(obj,t)
            %compute the approximated Hessian of the cost value
            [n_tp, n_state, n_contr] = getParams(obj);
            
            mesh        = obj.dyn.environment.mesh;
            
            alpha_   = obj.alpha;
            beta_    = obj.beta;
            gamma_   = obj.gamma;
            kappa_   = obj.kappa;
            
            if t >= n_tp
                CostD_t = [beta_*ones(3,1); ...
                           gamma_*ones(4,1); ...
                           kappa_*ones(6,1); ...
                           zeros(n_contr,1)];
            else
                CostD_t = [beta_*ones(3,1); ...
                           gamma_*ones(4,1); ...
                           kappa_*ones(6,1); ...
                           alpha_ * mesh(t)*ones(n_contr,1)];
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
            
            % Parametrisieren     
            o.dyn = cBQD;
            o.vec = o.dyn.vec;
            o.alpha = 2;
            o.beta = 3;
            o.gamma = 42;
            o.kappa = 1337;            
            o.timepoint = 23;
            o.cam_pos = @(timepoint) [ 2; 0; 5];
            
        end
    end
    
    
    
end

