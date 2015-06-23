classdef MultiShooting < TestEnv
    %   FORWEULER providing discretized ODE constraint using forward euler
    %   Use the multiple shooting approach to rewrite the ode constraint of the
    %   ocp as a number of equality constraints. This version here is based on
    %   one step of forward euler for each interval only.
    
    properties
        dyn;   % handle for the classDyn object providing the right hand side of the ODE
    end
    
    properties(Dependent)
        solver;
    end
    
    methods
        function res = get.solver(obj)
            res = obj.dyn.solver;
        end
        
    end
    
    
    
    methods
        %constructor
        function mS = MultiShooting(varargin)
            if(nargin == 0)
                global TEST;
                
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
                % constructor based on two input values
                % a classDyn element and a classOCPparam element
            elseif (nargin == 1)
                if (isa(varargin{1},'Dyn'))
                    mS.dyn = varargin{1};
                else
                    error('wrong class type for robot dynamics');
                end
            else
                error('wrong number of inputs');
            end
        end
        
        % other functions
        function [H, HD] = h(obj)
            % compute the equality constraints using forward euler
            [n_int, n_state, n_contr, mesh, n_var] = obj.getParams();
                                    
            if ((size(obj.dyn.contr,2)==n_int+1) && (size(obj.dyn.state,2)==n_int+1)...
                    &&(n_state == size(obj.dyn.state,1)) ...
                    &&(n_contr == size(obj.dyn.contr,1)))
                
                H        = zeros(n_int*n_state,1);
                state_val   = obj.dyn.state;
                
                % use vector notation to generate sparse matrix
                rvec = ones(1,n_int*(2*n_state+n_state*(n_state+n_contr)));
                cvec = ones(1,n_int*(2*n_state+n_state*(n_state+n_contr)));
                vvec = zeros(1,n_int*(2*n_state+n_state*(n_state+n_contr)));
                
                ind = 0;
                [old_intervals] = obj.solver.preToDo();
                for timepoint=1:n_int
                    [F, J] = obj.solver.ode(timepoint);
                    % Bestimme h
                    H((timepoint-1)*n_state+1:timepoint*n_state) = F - state_val(:,timepoint+1);
                    
                    % Bestimme hD
                    
%                     rvec(ind+1:ind+n_state)   = (timepoint-1)*n_state+1:timepoint*n_state;
%                     cvec(ind+1:ind+n_state)   = (timepoint-1)*(n_var)+1:...
%                         (timepoint-1)*(n_var)+n_state;
%                     vvec(ind+1:ind+n_state)   = ones(1,n_state);
%                     ind                     = ind + n_state;
                    
                     
                    [si,sj,sv]              = find(J);
                    sn                      = nnz(J);
                     
                    rvec(ind+1:ind+sn)      = (timepoint-1)*n_state+si;
                    cvec(ind+1:ind+sn)      = (timepoint-1)*(n_state+n_contr) + sj;
                    vvec(ind+1:ind+sn)      = sv;
                    ind                     = ind + sn;
                     
                    rvec(ind+1:ind+n_state)   = (timepoint-1)*n_state+1:timepoint*n_state;
                    cvec(ind+1:ind+n_state)   = timepoint*(n_var)+1:...
                                     timepoint*(n_var)+n_state;
                    vvec(ind+1:ind+n_state)   = -ones(1,n_state);
                    ind                     = ind + n_state;
                end
                HD = sparse(rvec(1:ind),cvec(1:ind),vvec(1:ind),...
                      n_int*n_state,(n_int+1)*(n_var));
                obj.solver.postToDo(old_intervals);
            else
                error('wrong state and control lengths wrt index.');
            end
        end
        
%         %function val = hD(obj)
%             % compute Jacobian of the equality constraints using forward euler
%             [n_int, n_state, n_contr, mesh] = getParams(obj);
%             
%             
%             if ((size(obj.dyn.contr,2)==n_int+1) && (size(obj.dyn.state,2)==n_int+1)...
%                     &&(n_state == size(obj.dyn.state,1)) ...
%                     &&(n_contr == size(obj.dyn.contr,1)))
%                 
%                 error('Inkonstruktion')
%                 
% %                 %val = sparse(n_int*n_var,(n_int+1)*(n_var+n_contr));
% %                 % use vector notation to generate sparse matrix
% %                 rvec = ones(1,n_int*(2*n_state+n_state*(n_state+n_contr)));
% %                 cvec = ones(1,n_int*(2*n_state+n_state*(n_state+n_contr)));
% %                 vvec = zeros(1,n_int*(2*n_state+n_state*(n_state+n_contr)));
% %                 
% %                 ind = 0;
% %                 
% %                 for i=1:n_int
% %                     rvec(ind+1:ind+n_state)   = (i-1)*n_state+1:i*n_state;
% %                     cvec(ind+1:ind+n_state)   = (i-1)*(n_state+n_contr)+1:...
% %                         (i-1)*(n_state+n_contr)+n_state;
% %                     vvec(ind+1:ind+n_state)   = ones(1,n_state);
% %                     ind                     = ind + n_state;
% %                     
% %                     mat                     = mesh(i)*obj.dyn.dotD(i);
% %                     [si,sj,sv]              = find(mat);
% %                     sn                      = length(sv);
% %                     
% %                     rvec(ind+1:ind+sn)      = (i-1)*n_state+si;
% %                     cvec(ind+1:ind+sn)      = (i-1)*(n_state+n_contr) + sj;
% %                     vvec(ind+1:ind+sn)      = sv;
% %                     ind                     = ind + sn;
% %                     
% %                     rvec(ind+1:ind+n_state)   = (i-1)*n_state+1:i*n_state;
% %                     cvec(ind+1:ind+n_state)   = i*(n_state+n_contr)+1:...
% %                         i*(n_state+n_contr)+n_state;
% %                     vvec(ind+1:ind+n_state)   = -ones(1,n_state);
% %                     ind                     = ind + n_state;
% %                 end
% %                 
% %                 val = sparse(rvec(1:ind),cvec(1:ind),vvec(1:ind),...
% %                     n_int*n_state,(n_int+1)*(n_state+n_contr));
%             else
%                 error('wrong state and control lengths wrt index.');
%             end
%         end
        %Approximieren durch Euler
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
        
        function [n_int, n_state, n_contr, mesh, n_var] = getParams(obj)
            n_int       = obj.dyn.environment.n_intervals;
            n_state     = obj.dyn.robot.n_state;
            n_contr     = obj.dyn.robot.n_contr;
            mesh        = obj.dyn.environment.mesh;
            n_var     = obj.dyn.robot.n_var;
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
            FE = ForwEuler();
            
            cBQD = BasisQDyn(model, env, FE);
            cBQD.vec = rand(model.n_var * (n_intervals+1),1);
            obj.dyn = cBQD;
        end
        
        function [vec_old, n, m, n_timepoints, dyn] = setup(obj, func)
            vec_old = obj.dyn.vec;
            n_timepoints = obj.dyn.environment.n_timepoints;
            dyn = obj.dyn;
            n = obj.dyn.robot.n_var;
            m = size(func());
            m = m(1);
        end
        
        function hD =  gethD(obj)
            [h, hD] = obj.h();
        end
        
        
    end
    
    methods(Test)
        function testhD(obj)
            % TESTHD This method derives numerically obj.h and compares it
            %with obj.hD
            n_intervals = 2;
            obj.setupTest(n_intervals);
            [n_int, n_state, n_contr, mesh] = getParams(obj);
            
            [h, anaDiff] = obj.h();
            func = @() obj.h;
            numDiff = obj.numDiff_nD_AllT(func);
            
            obj.assertSize(anaDiff, size(numDiff) );
            obj.assertSize(anaDiff, [n_intervals * 13 , (n_intervals+1)* 17 ]);
            obj.assertLessThan(max(abs(anaDiff - numDiff)), obj.tol);
        end
        
        function testhDD(obj)
            % TESTHDD This methods derives numerically obj.hd and compares
            % it with obj.hDD
            n_intervals = 3;
            obj.setupTest(n_intervals);
            [n_int, n_state, n_contr, mesh, n_var] = getParams(obj);
            
            func = @() obj.gethD;
            numDiff = obj.numDiff_nxnD_AllT(func);
            anaDiff = obj.hDD();
            
            size_nDiff_i = (n_var) * (n_intervals +1 );
            for i = 1:(n_int * n_state)
                numDiff_i = reshape(numDiff(i,:,:), [size_nDiff_i size_nDiff_i]);
                obj.assertSize(anaDiff{1}, size(numDiff_i));
                obj.assertLessThan(max(abs(anaDiff{i} - numDiff_i)), obj.tol);
            end
        end
    end
end