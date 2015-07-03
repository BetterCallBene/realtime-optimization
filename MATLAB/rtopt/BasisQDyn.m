classdef BasisQDyn < BasisGenQDyn
    
    properties
        backdoor_vec;
        
        flagSteadyPoint;
        
        dotDpatternflag;
        dotDDpatternflag;
        
        dotDCellpattern;
        dotDDCellpattern;
        
        steadyPoint;
    end
    
    properties(Dependent)
        state;
        contr;
        vec;
    end
    
    
    methods
        function bQDyn = BasisQDyn(varargin)
            bQDyn@BasisGenQDyn();
            
            bQDyn.dotDpatternflag = true;
            bQDyn.dotDDpatternflag = true;
            bQDyn.flagSteadyPoint = true;
            
            nargin_tmp = nargin;
            m = metaclass(bQDyn);
            
            if strcmp(m.Name, 'BasisQDyn') == false %Supclass?
                varargin = varargin{1};
                nargin_tmp = length(varargin);
            end
            
            if (nargin_tmp >= 1)
                mc = metaclass(varargin{1});
                if strcmp(mc.SuperclassList(1).Name, 'Model')
                    bQDyn.robot = varargin{1};
                else
                    error('First argument must be instance of Model');
                end
            end
            if (nargin_tmp >= 2)
                if isa(varargin{2}, 'Environment')
                    bQDyn.environment = varargin{2};
                else
                    error('Second argument must be instance of Environment');
                end
            end
            
            if (nargin_tmp >= 3)
                if isa(varargin{3}, 'Solver')
                    bQDyn.solver = varargin{3};
                    bQDyn.solver.dyn = bQDyn;
                else
                    error('Third argument must be instance of Solver');
                end
            end
        end
        
        %set methods
        function set.backdoor_vec(obj,vec)
            if (isempty(obj.backdoor_vec) || size(obj.backdoor_vec, 1) ~= size(vec, 1) || (norm(obj.backdoor_vec - vec)>1e-10))
                obj.backdoor_vec = vec;
                obj.emptyResults();
            end
        end
        
        function res = get.backdoor_vec(obj)
            res = obj.backdoor_vec;
        end
        
        function set.vec(obj, vec)
            n_int = obj.environment.n_intervals;
            n_state = obj.robot.n_state;
            n_contr = obj.robot.n_contr;
            if isempty(n_int)
                error('Gitter ist noch nicht initialisiert');
            end
            n_var = (n_state + n_contr);
            n = (n_int + 1)*n_var;
            if size(vec, 1) == n
                for i = 1:(n_int+1)
                    q = vec((i-1)*n_var + 4:(i-1)*n_var + 7);
                    vec((i-1)*n_var + 4:(i-1)*n_var + 7) = q./norm(q);
                end
                obj.backdoor_vec = vec;
            else
                error(strcat('Gr��e der State Matrix ist falsch. Erwartete Gr��e  ',int2str((n_int + 1)*(n_state + n_contr)),'xn'));
            end
        end
        
        function ret =  get.vec(obj)
            ret = obj.backdoor_vec;
        end
        
        % getter methods
        function val = get.state(obj)
            % get the current state values
            % (if not yet stored, extract them from vec)
            
            n_int       = obj.environment.n_intervals;
            n_state       = obj.robot.n_state;
            n_contr     = obj.robot.n_contr;
            
            val        = zeros(n_state,n_int+1);
            
            for i = 1:n_int+1
                val(:,i) = obj.vec((i-1)*(n_state+n_contr)+1:...
                    (i-1)*(n_state+n_contr)+n_state);
            end
            
        end
        
        
        function val = get.contr(obj)
            % get the current control values
            % (if not yet stored, extract them from vec)
            
            n_int       = obj.environment.n_intervals;
            n_var       = obj.robot.n_state;
            n_contr     = obj.robot.n_contr;
            
            val        = zeros(n_contr,n_int+1);
            
            for i = 1:n_int+1
                val(:,i) = obj.vec((i-1)*(n_var+n_contr)+n_var + 1:...
                    i*(n_var+n_contr));
            end
        end
        
        function ret = get.steadyPoint(obj)
            if obj.flagSteadyPoint
                %[q,v,omega,u,Iges,IM,m,kT,kQ,d,g] = obj.getParams();
                m = obj.robot.m;
                kT = obj.robot.kT;
                g = obj.environment.g;
                estimator_u = sqrt(1/4 * m * g * 1/kT);
                estimator = [zeros(3, 1);1;zeros(9, 1);repmat(estimator_u, 4, 1)];
                options = optimoptions('fsolve', 'Algorithm', 'levenberg-marquardt');
                n_int_sav = obj.environment.n_intervals;
                obj.environment.n_intervals = 0;
                [obj.steadyPoint,fval,exitflag,output] = fsolve(@obj.helperF, estimator, options);
                obj.environment.n_intervals = n_int_sav ;
                obj.flagSteadyPoint = false;
            end
            ret = obj.steadyPoint;
        end
        
        function res = helperF(obj, y)
            obj.backdoor_vec = y;
            res = obj.F(:, 1);
        end
        
        
        function res = dot(obj,ind) % Rechte Seite der ODE aus (Quadrocoptermodell.pdf: (1.46))
            % compute the right hand side of the ode for a given time
            % instance ind
            res = obj.F(:, ind);
        end
        
        
        function res = dotD(obj,ind )
            % compute the Jacobian of the right hand side of the ode for
            % a given time instance ind
            J_ = obj.J;
            res = sparse(J_(:, 1), J_(:, 2), J_(:, ind + 2), 13, 17);
        end
        
        %         function res = dotDD(obj,ind)
        %             % compute the Hessian of the right hand side of the ode for
        %             % a given time instance ind
        %             n_state = obj.robot.n_state;
        %             res = cell(1, n_state);
        %             H_ = obj.H;
        %             for i = 1:size(H_, 1)
        %                 if isempty(res{H_(i, 1)})
        %                     res{H_(i, 1)} = sparse(17, 17);
        %                 end
        %                 res{H_(i, 1)}(H_(i, 2), H_(i, 3)) = H_(i, ind + 3);
        %             end
        %         end
        
        function res = dotDD(obj,ind)
            % compute the Hessian of the right hand side of the ode for
            %           % a given time instance ind
            res = zeros(13, 17, 17);
            H_ = obj.H;
            for i = 1:size(H_, 1)
                res(H_(i, 1), H_(i, 2), H_(i, 3)) = H_(i, ind + 3);
            end
        end
    end
    methods
        function pattern =  dotDpattern(obj)
            if obj.dotDpatternflag
                J_ = obj.J;
                
                obj.dotDCellpattern = sparse(J_(:, 1), J_(:, 2), ones(size(J_(:, 1), 1), size(J_(:, 1), 2)), 13, 17);
                obj.dotDpatternflag = false;
            end
            pattern = obj.dotDCellpattern;
        end
        
        function pattern =  dotDDpattern(obj)
            if obj.dotDDpatternflag
                n_state       = obj.robot.n_state;
                
                obj.dotDDCellpattern = cell(1, n_state);
                H_ = obj.H;
                for i = 1:size(H_, 1)
                    if isempty(obj.dotDDCellpattern{H_(i, 1)})
                        obj.dotDDCellpattern{H_(i, 1)} = sparse(17, 17);
                    end
                    obj.dotDDCellpattern{H_(i, 1)}(H_(i, 2), H_(i, 3)) = 1;
                end
                
                obj.dotDDpatternflag = false;
            end
            pattern = obj.dotDDCellpattern;
        end
        
    end
    methods(Test)
        function testInitBasisQDyn(obj)
            n_int_ = uint16(50);
            obj.environment = Environment();
            obj.environment.setUniformMesh(n_int_);
            
            obj.robot = Quadrocopter();
            n_state_ = obj.robot.n_state;
            n_contr_ = obj.robot.n_contr;
            
            n_var = n_state_ + n_contr_;
            obj.vec = rand(n_var * (n_int_ + 1), 1);
            %obj.state = rand(n_state_, n_int_+1);
            %obj.contr = rand(n_contr_, n_int_+1);
            
            %for ind=1:n_int_
            %    var1 = obj.dot(ind);
            %    var2 = obj.dotD(ind);
            %    var3 = obj.dotDD(ind);
            %end
        end
    end
end