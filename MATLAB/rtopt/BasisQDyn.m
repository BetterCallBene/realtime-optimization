classdef BasisQDyn < BasisGenQDyn
    
    properties
        backdoor_vec;
    end
    
    properties(Dependent)
        state;
        contr;
        vec;
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
        
        %set methods
        function set.backdoor_vec(obj,vec)
            % set new input vector 
            % and pass information on th dode 
            if (isempty(obj.backdoor_vec) || (norm(obj.backdoor_vec - vec)>1e-10))
                obj.backdoor_vec = vec;
                obj.emptyResults();
            end
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
        
        %function res = get.state(cq)
        %    res = cq.backdoor_state;
        %end
        
        
        
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
                for i = 1:size(J_, 1)
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
                for i = 1:size(H_, 1)
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