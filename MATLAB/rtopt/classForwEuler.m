classdef classForwEuler < handle
% classForwEuler providing discretized ODE constraint using forward euler
%   Use the multiple shooting approach to rewrite the ode constraint of the
%   ocp as a number of equality constraints. This version here is based on
%   one step of forward euler for each interval only. 

    properties
        state; % the state vector combining the joint angles theta and their first time derivates dot_theta for all time instances          
        contr; % the control vector combining the two joint torque values for all time instances
        dyn;   % handle for the classDyn object providing the right hand side of the ODE
        param; % handel for the classOCPparam object providing the properties describing the ocp
    end
    
    methods
        %constructor
        function cFE = classForwEuler(varargin)
            % constructor based on two input values
            % a classDyn element and a classOCPparam element
            if (nargin == 2)
                if (isa(varargin{1},'classDyn'))
                    cFE.dyn = varargin{1};
                else
                    error('wrong class type for robot dynamics');
                end
                
                if (isa(varargin{2},'classOCPparam'))
                    cFE.param = varargin{2};
                else
                    error('wrong class type for problem parameter');
                end
            else
                error('wrong number of inputs');     
            end
        end
        
        %set methods
        function set.state(obj,state)
            % set the current values for the states at all time instances
            % and pass it on the dyn         
            obj.state = state;
            if (~isempty(obj.dyn))
                obj.dyn.state = state;
            end
        end
        
        function set.contr(obj,contr)
            % set the current values for the controls at all time instances
            % and pass it on the dyn                
            obj.contr = contr;
            if (~isempty(obj.dyn))
                obj.dyn.contr = contr;
            end
        end
        
        % other functions
        function val = h(obj)
            % compute the equality constraints using forward euler
            n_int       = obj.param.n_intervals;
            n_var       = obj.param.n_var;
            n_contr     = obj.param.n_contr;
            mesh        = obj.param.mesh;
            
            if ((size(obj.contr,2)==n_int+1) && (size(obj.state,2)==n_int+1)...
                    &&(n_var == size(obj.state,1)) ...
                    &&(n_contr == size(obj.contr,1)))

                    val         = zeros(n_int*n_var,1);
                    state_val   = obj.state;
                    
                    for i=1:n_int
                        val((i-1)*n_var+1:i*n_var) = state_val(:,i) + ...
                            mesh(i)*obj.dyn.dot(i)-state_val(:,i+1);
                    end
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function val = hD(obj)
             % compute Jacobian of the equality constraints using forward euler
            n_int       = obj.param.n_intervals;
            n_var       = obj.param.n_var;
            n_contr     = obj.param.n_contr;
            mesh        = obj.param.mesh;
            
            if ((size(obj.contr,2)==n_int+1) && (size(obj.state,2)==n_int+1)...
                    &&(n_var == size(obj.state,1)) ...
                    &&(n_contr == size(obj.contr,1)))

                    %val = sparse(n_int*n_var,(n_int+1)*(n_var+n_contr));
                    % use vector notation to generate sparse matrix
                    rvec = ones(1,n_int*(2*n_var+n_var*(n_var+n_contr)));
                    cvec = ones(1,n_int*(2*n_var+n_var*(n_var+n_contr)));
                    vvec = zeros(1,n_int*(2*n_var+n_var*(n_var+n_contr)));
                    
                    ind = 0;
                    
                    for i=1:n_int
                        rvec(ind+1:ind+n_var)   = (i-1)*n_var+1:i*n_var;
                        cvec(ind+1:ind+n_var)   = (i-1)*(n_var+n_contr)+1:...
                                    (i-1)*(n_var+n_contr)+n_var;
                        vvec(ind+1:ind+n_var)   = ones(1,n_var);
                        ind                     = ind + n_var;
                        
                        mat                     = mesh(i)*obj.dyn.dotD(i);
                        [si,sj,sv]              = find(mat);
                        sn                      = length(sv);
                        
                        rvec(ind+1:ind+sn)      = (i-1)*n_var+si;
                        cvec(ind+1:ind+sn)      = (i-1)*(n_var+n_contr) + sj;
                        vvec(ind+1:ind+sn)      = sv;
                        ind                     = ind + sn;
                        
                        rvec(ind+1:ind+n_var)   = (i-1)*n_var+1:i*n_var;
                        cvec(ind+1:ind+n_var)   = i*(n_var+n_contr)+1:...
                                    i*(n_var+n_contr)+n_var;
                        vvec(ind+1:ind+n_var)   = -ones(1,n_var);
                        ind                     = ind + n_var;                        
                    end
                    
                    val = sparse(rvec(1:ind),cvec(1:ind),vvec(1:ind),...
                                    n_int*n_var,(n_int+1)*(n_var+n_contr));
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function val = hDD(obj)
             % compute the Hessian the equality constraints using forward euler
            n_int       = obj.param.n_intervals;
            n_var       = obj.param.n_var;
            n_contr     = obj.param.n_contr;
            mesh        = obj.param.mesh;
            
            if ((size(obj.contr,2)==n_int+1) && (size(obj.state,2)==n_int+1)...
                    &&(n_var == size(obj.state,1)) ...
                    &&(n_contr == size(obj.contr,1)))
                
                val     = cell(n_int*n_var,1);

                for i=1:n_int
                    dotDD       = obj.dyn.dotDD(i);
                    
                    for j=1:length(dotDD)
                        [si,sj,sv]  = find(mesh(i)*dotDD{j});                       
                        si          = si + (i-1)*(n_var+n_contr);
                        sj          = sj + (i-1)*(n_var+n_contr);

                        val{(i-1)*n_var+j} = sparse(si,sj,sv,...
                            (n_int+1)*(n_var+n_contr),(n_int+1)*(n_var+n_contr));
                    end
                end
            else
                error('wrong state and control lengths wrt index.');
            end
        end       
        
    end
end