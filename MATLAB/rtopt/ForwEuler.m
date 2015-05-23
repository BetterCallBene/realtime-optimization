classdef ForwEuler < TestEnv
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
            
            if(nargin == 0)
                global TEST;
                
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
                
                
                
                
                % constructor based on two input values
                % a classDyn element and a classOCPparam element
            elseif (nargin == 1)
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
            
            if ((size(obj.dyn.contr,2)==n_int+1) && (size(obj.dyn.state,2)==n_int+1)...
                    &&(n_state == size(obj.dyn.state,1)) ...
                    &&(n_contr == size(obj.dyn.contr,1)))
                
                val         = zeros(n_int*n_state,1);
                state_val   = obj.dyn.state;
                
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
            
            
            if ((size(obj.dyn.contr,2)==n_int+1) && (size(obj.dyn.state,2)==n_int+1)...
                    &&(n_state == size(obj.dyn.state,1)) ...
                    &&(n_contr == size(obj.dyn.contr,1)))
                
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
            [n_int, n_state, n_contr, mesh] = getParams(obj);
            
            if ((size(obj.dyn.contr,2)==n_int+1) && (size(obj.dyn.state,2)==n_int+1)...
                    &&(n_state == size(obj.dyn.state,1)) ...
                    &&(n_contr == size(obj.dyn.contr,1)))
                
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
            n_state     = obj.dyn.robot.n_state;
            n_contr     = obj.dyn.robot.n_contr;
            mesh        = obj.dyn.environment.mesh;
        end
        
        function setupTest(obj,n_intervals)
            % Quadrocopter soll 5 Meter hoch fliegen
            xbc = [         ... Variablenname L�nge   Name
                ... Anfangsbedingung
                0, 0, 0,    ...     r           3      Ortsvektor
                1, 0, 0, 0, ...     q           4      Quaternion (Einheitsquaternion)
                0, 0, 0,    ...     v           3      Translatorische Geschwindigkeit
                0, 0, 0;    ...     w           3      Winkelgeschwindigkeit
                ... Endbedingung
                0, 0, 5,    ...
                1, 0, 0, 0, ...
                0, 0, 0,    ...
                0, 0, 0     ...
                ];
            
            env = Environment();
            env.xbc = xbc;
            env.setUniformMesh(uint8(n_intervals));
            
            model = Quadrocopter();
            
            cBQD = BasisQDyn(model, env);
            cBQD.vec = rand(model.n_var * (n_intervals+1),1);
            obj.dyn = cBQD;
            
        end
        
        function [vec_old, n,m,n_timepoints] = setup(obj,func)
            vec_old = obj.dyn.vec;
            n_timepoints = obj.dyn.environment.n_timepoints;
            n = obj.dyn.robot.n_var;
            m = size(func());
            m=m(1);
        end
        
        function func_p = plusEpsShift(obj, i ,t ,vec_old,func)
            vec_p = vec_old;
            vec_p((t-1) * obj.dyn.robot.n_var + i ) = vec_p((t-1) * obj.dyn.robot.n_var + i) + obj.eps;
            obj.dyn.backdoor_vec = vec_p;
            func_p = func();
            obj.dyn.vec = vec_old;
        end
        
        function func_n = minusEpsShift(obj, i ,t ,vec_old,func)
            vec_n = vec_old;
            vec_n((t-1) * obj.dyn.robot.n_var + i ) = vec_n((t-1) * obj.dyn.robot.n_var + i) - obj.eps;
            obj.dyn.backdoor_vec = vec_n;
            func_n = func();
            obj.dyn.vec = vec_old;
        end
    end
    
    methods(Test)
        function testhD(obj)
            % TESTHD This method derives numerically obj.h and compares it
            %with obj.hD
            n_intervals = 25;
            obj.setupTest(n_intervals);
            [n_int, n_state, n_contr, mesh] = getParams(obj);
            
            func = @() obj.h;
            numDiff = obj.numDiff_nD_AllT(func);
            anaDiff = obj.hD();
            
            obj.assertSize(anaDiff, size(numDiff) );
            obj.assertSize(anaDiff, [n_intervals * 13 , (n_intervals+1)* 17 ]);
            obj.assertLessThan(max(abs(anaDiff - numDiff)), obj.tol);
            
        end
        
        function testhDD(obj)
            % TESTHDD This methods derives numerically obj.hd and compares
            % it with obj.hDD
            n_intervals = 25;
            obj.setupTest(n_intervals);
            [n_int, n_state, n_contr, mesh] = getParams(obj);
            
            func = @() obj.hD;
            numDiff = obj.numDiff_nxnD_AllT(func);
            anaDiff = obj.hDD();
            
            size_nDiff_i = (n_state + n_contr) * (n_intervals +1 );
            for i = 1:(n_intervals * n_state)
                numDiff_i = reshape(numDiff(i,:,:), [size_nDiff_i size_nDiff_i]);
                obj.assertSize(anaDiff{1}, size(numDiff_i));
                obj.assertLessThan(max(abs(anaDiff{i} - numDiff_i)), obj.tol);
            end
        end
    end
end