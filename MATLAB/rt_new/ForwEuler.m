classdef ForwEuler < handle
    %   FORWEULER providing discretized ODE constraint using forward euler
    %   Use the multiple shooting approach to rewrite the ode constraint of the
    %   ocp as a number of equality constraints. This version here is based on
    %   one step of forward euler for each interval only.
    
    properties
        dyn;   % handle for the classDyn object providing the right hand side of the ODE
    end
    
    methods
        %constructor
        function fE = ForwEuler(varargin)
            % constructor based on two input values
            % a classDyn element and a classOCPparam element
            if (nargin == 1)
                if (isa(varargin{1},'Dyn'))
                    fE.dyn = varargin{1};
                else
                    error('wrong class type for robot dynamics');
                end
                
            else
                error('wrong number of inputs');
            end
        end
        
        % other functions
        function val = h(obj)
            % compute the equality constraints using forward euler
            n_int       = obj.dyn.environment.n_intervals;
            n_state       = obj.dyn.robot.n_state;
            n_contr     = obj.dyn.robot.n_contr;
            mesh        = obj.dyn.environment.mesh;
            
            if ((size(obj.contr,2)==n_int+1) && (size(obj.state,2)==n_int+1)...
                    &&(n_state == size(obj.state,1)) ...
                    &&(n_contr == size(obj.contr,1)))
                
                val         = zeros(n_int*n_state,1);
                state_val   = obj.state;
                
                for i=1:n_int
                    val((i-1)*n_state+1:i*n_state) = state_val(:,i) + ...
                        mesh(i)*obj.dyn.dot(i)-state_val(:,i+1);
                end
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function val = hD(obj)
            % compute Jacobian of the equality constraints using forward euler
            [n_int, n_state, n_contr, mesh] = getParams(obj);
            
            
            if ((size(obj.contr,2)==n_int+1) && (size(obj.state,2)==n_int+1)...
                    &&(n_state == size(obj.state,1)) ...
                    &&(n_contr == size(obj.contr,1)))
                
                %val = sparse(n_int*n_var,(n_int+1)*(n_var+n_contr));
                % use vector notation to generate sparse matrix
                rvec = ones(1,n_int*(2*n_state+n_state*(n_state+n_contr)));
                cvec = ones(1,n_int*(2*n_state+n_state*(n_state+n_contr)));
                vvec = zeros(1,n_int*(2*n_state+n_state*(n_state+n_contr)));
                
                ind = 0;
                
                for i=1:n_int
                    rvec(ind+1:ind+n_state)   = (i-1)*n_state+1:i*n_state;
                    cvec(ind+1:ind+n_state)   = (i-1)*(n_state+n_contr)+1:...
                        (i-1)*(n_state+n_contr)+n_state;
                    vvec(ind+1:ind+n_state)   = ones(1,n_state);
                    ind                     = ind + n_state;
                    
                    mat                     = mesh(i)*obj.dyn.dotD(i);
                    [si,sj,sv]              = find(mat);
                    sn                      = length(sv);
                    
                    rvec(ind+1:ind+sn)      = (i-1)*n_state+si;
                    cvec(ind+1:ind+sn)      = (i-1)*(n_state+n_contr) + sj;
                    vvec(ind+1:ind+sn)      = sv;
                    ind                     = ind + sn;
                    
                    rvec(ind+1:ind+n_state)   = (i-1)*n_state+1:i*n_state;
                    cvec(ind+1:ind+n_state)   = i*(n_state+n_contr)+1:...
                        i*(n_state+n_contr)+n_state;
                    vvec(ind+1:ind+n_state)   = -ones(1,n_state);
                    ind                     = ind + n_state;
                end
                
                val = sparse(rvec(1:ind),cvec(1:ind),vvec(1:ind),...
                    n_int*n_state,(n_int+1)*(n_state+n_contr));
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function val = hDD(obj)
            % compute the Hessian the equality constraints using forward euler
            [n_int, n_state, n_contr, mesh] = getParams(obj)
            
            if ((size(obj.contr,2)==n_int+1) && (size(obj.state,2)==n_int+1)...
                    &&(n_state == size(obj.state,1)) ...
                    &&(n_contr == size(obj.contr,1)))
                
                val     = cell(n_int*n_state,1);
                
                for i=1:n_int
                    dotDD       = obj.dyn.dotDD(i);
                    
                    for j=1:length(dotDD)
                        [si,sj,sv]  = find(mesh(i)*dotDD{j});
                        si          = si + (i-1)*(n_state+n_contr);
                        sj          = sj + (i-1)*(n_state+n_contr);
                        
                        val{(i-1)*n_state+j} = sparse(si,sj,sv,...
                            (n_int+1)*(n_state+n_contr),(n_int+1)*(n_state+n_contr));
                    end
                end
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
        function [n_int, n_state, n_contr, mesh] = getParams(obj)
            n_int       = obj.dyn.environment.n_intervals;
            n_state       = obj.dyn.robot.n_state;
            n_contr     = obj.dyn.robot.n_contr;
            mesh        = obj.dyn.environment.mesh;
        end
        
    end
end