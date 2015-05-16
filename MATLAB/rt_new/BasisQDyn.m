classdef BasisQDyn < BasisGenQDyn
    
    properties
        backdoor_state; % the state matrix [r, q, v, w] in R^(13xn) for all time instances
        contr; % the control matrix [u1, u2, u3, u4] in R^(4xn) for all time instances
    end
    
    properties(Dependent)
        state;
    end
    
    methods
        function bQDyn = BasisQDyn(varargin)
            bQDyn@BasisGenQDyn();
            if (nargin >= 1)
                mc = metaclass(varargin{1});
                if strcmp(mc.SuperclassList(1).Name, 'Model')
                    bQDyn.robot = varargin{1};
                else
                    error('First argument must be instance of Model');
                end
            end
            if (nargin >= 2)
                if isa(varargin{2}, 'Environment')
                    bQDyn.environment = varargin{2};
                else
                    error('Second argument must be instance of Environment');
                end
            end
        end
        
        function set.state(cq, state)
            
            if isempty(cq.environment.n_timepoints)
                error('Gitter ist noch nicht initialisiert');
            end
            
            %   if size(state, 1) == cq.robot.n_state && size(state, 2) == cq.environment.n_timepoints;
            
            n = size(state, 2);
            for i=1:n
                q = state(4:7, i);
                state(4:7, i) = q./norm(q);
            end
            
            cq.backdoor_state = state;
            cq.emptyResults();
            %  else
            %                 error(strcat('Gr��e der State Matrix ist falsch. Erwartete Gr��e  ',int2str(cq.n_state),'xn'));
            % end
        end
        
        
        function set.backdoor_state(cq, state)
            cq.backdoor_state = state;
        end
        
        function res = get.state(cq)
            res = cq.backdoor_state;
        end
        
        function set.contr(cq, cntrl)
            
            if isempty(cq.environment.n_timepoints)
                error('Gitter ist noch nicht initialisiert');
            end
            
            %   if size(cntrl, 1) == cq.robot.n_contr && size(cntrl, 2) == cq.environment.n_timepoints
            
            cq.contr = cntrl;
            cq.emptyResults();
            
            %else
            %     error(strcat('Gr��e der State Matrix ist falsch. Erwartete Gr��e  ',int2str(cq.n_contr),'xn'));
            %end
        end
        
        function res = dot(obj,ind) % Rechte Seite der ODE aus (Quadrocoptermodell.pdf: (1.46))
            % compute the right hand side of the ode for a given time
            % instance ind
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                res = obj.F(:, ind);
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        
        function res = dotD(obj,ind)
            % compute the Jacobian of the right hand side of the ode for
            % a given time instance ind
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                
                J_ = obj.J;
                res = zeros(13, 17);
                for i = 1:length(J_)
                    res(J_(i, 1), J_(i, 2)) = J_(i, ind + 2);
                end
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function res = dotDD(obj,ind)
            % compute the Hessian of the right hand side of the ode for
            % a given time instance ind
            if ((size(obj.contr,2)==size(obj.state,2))&&(ind <= size(obj.contr,2)))
                res = cell(1, 13);
                H_ = obj.H;
                for i = 1:length(H_)
                    if isempty(res{H_(i, 1)})
                        res{H_(i, 1)} = zeros(17, 17);
                    end
                    res{H_(i, 1)}(H_(i, 2), H_(i, 3)) = H_(i, ind + 3);
                end
            end
        end
        
    end
    methods(Test)
        function testInitBasisQDyn(obj)
            n_int_ = uint8(50);
            obj.environment = Environment();
            obj.environment.setUniformMesh(n_int_);
            
            obj.robot = Quadrocopter();
            n_state_ = obj.robot.n_state;
            n_contr_ = obj.robot.n_contr;
            
            obj.state = rand(n_state_, n_int_+1);
            obj.contr = rand(n_contr_, n_int_+1);
            
            for ind=1:n_int_
                var1 = obj.dot(ind);
                var2 = obj.dotD(ind);
                var3 = obj.dotDD(ind);
            end
        end
    end
end