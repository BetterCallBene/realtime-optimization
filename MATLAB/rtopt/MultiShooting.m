classdef MultiShooting < TestEnv
    %   FORWEULER providing discretized ODE constraint using forward euler
    %   Use the multiple shooting approach to rewrite the ode constraint of the
    %   ocp as a number of equality constraints. This version here is based on
    %   one step of forward euler for each interval only.
    
    properties
        dyn;   % handle for the classDyn object providing the right hand side of the ODE
        flag_h;
        flag_hDD;
        
        cacheH;
        cacheHD;
        cacheHDD;
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
            
            if obj.flag_h
                H = obj.cacheH;
                HD = obj.cacheHD;
            else
                
                
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
                    
                    Fs = cell(n_int, 1);
                    Js = cell(n_int, 1);
                    %Parallelauswertung der odes
                    parfor timepoint=1:n_int
                        [F, J] = obj.solver.ode(timepoint);
                        Fs{timepoint} = F;
                        Js{timepoint} = J;
                    end
                    %Sortierung der F, J
                    for timepoint=1:n_int
                        %[F, J] = obj.solver.ode(timepoint);
                        F = Fs{timepoint};
                        J = Js{timepoint};
                        % Bestimme h
                        H((timepoint-1)*n_state+1:timepoint*n_state) = F - state_val(:,timepoint+1);
                        
                        % Bestimme hD
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
                    
                    obj.cacheH = H;
                    obj.cacheHD = HD;
                    obj.flag_h  = true;
                    
                else
                    error('wrong state and control lengths wrt index.');
                    
                    Fs = cell(n_int, 1);
                    Js = cell(n_int, 1);
                    %Parallelauswertung der odes
                    parfor timepoint=1:n_int
                        [F, J] = obj.solver.ode(timepoint);
                        Fs{timepoint} = F;
                        Js{timepoint} = J;
                    end
                    %Sortierung der F, J
                    for timepoint=1:n_int
                        %[F, J] = obj.solver.ode(timepoint);
                        F = Fs{timepoint};
                        J = Js{timepoint};
                        % Bestimme h
                        
                        H((timepoint-1)*n_state+1:timepoint*n_state) = F - state_val(:,timepoint+1);
                        
                        % Bestimme hD
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
                end
            end
            
            function HDD = hDD(obj)
                if obj.flag_hDD
                    HDD = obj.cacheHDD;
                else
                    % compute the Hessian the equality constraints using forward euler
                    [n_int, n_state, n_contr, mesh] = getParams(obj);
                    
                    disp('for testing only');
                    
                    if ((size(obj.dyn.contr,2)==n_int+1) && (size(obj.dyn.state,2)==n_int+1)...
                            &&(n_state == size(obj.dyn.state,1)) ...
                            &&(n_contr == size(obj.dyn.contr,1)))
                        
                        HDD     = cell(n_int*n_state,1);
                        
                        for i=1:n_int
                            dotDD       = obj.dyn.dotDD(i);
                            
                            for j=1:length(dotDD)
                                [si,sj,sv]  = find(mesh(i)*dotDD{j});
                                si          = si + (i-1)*(n_state+n_contr);
                                sj          = sj + (i-1)*(n_state+n_contr);
                                
                                HDD{(i-1)*n_state+j} = sparse(si,sj,sv,...
                                    (n_int+1)*(n_state+n_contr),(n_int+1)*(n_state+n_contr));
                            end
                        end
                        
                        obj.cacheHDD = HDD;
                        obj.flag_hDD = true;
                    else
                        error('wrong state and control lengths wrt index.');
                    end
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
                opts_ = odeset('RelTol',1e-4,'AbsTol',1e-5);
                integrator = ode15sM(opts_);
                
                %load('TestData', 'data');
                cBQD = BasisQDyn(model, env, integrator);
                cBQD.vec = rand(17 * (n_intervals+1), 1);
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
            
            function numDiff = numDiff_nD_AllT(obj, func)
                % NUMDIFF_NDALLT This method calculates numericalle the
                % derivative of func, but for all timepoints
                
                [vec_old, n, m, n_timepoints, dyn] = obj.setup(func);
                numDiff = zeros(m , n * n_timepoints);
                
                for timepoint = 1:n_timepoints
                    for i = 1:n
                        func_p = obj.plusEpsShift(i,timepoint,vec_old, func, n, dyn);
                        func_n = obj.minusEpsShift(i,timepoint,vec_old, func, n, dyn);
                        
                        %Central difference
                        numDiff( : , ((timepoint-1) * n) + i) ...
                            = (func_p - func_n)/2/obj.eps;
                    end
                end
            end
            
            function hD =  gethD(obj)
                [h, hD] = obj.h();
            end
            
            
            
            
        end
        
        methods(Test)
            
            function testh(testCase)
                n_intervals = 20;
                testCase.setupTest(n_intervals);
                [n_int, n_state, n_contr, mesh] = getParams(testCase);
                tic
                [h, anaDiff] = testCase.h();
                toc
            end
            
            function testhD(testCase)
                n_intervals = 3;
                testCase.setupTest(n_intervals);
                [n_int, n_state, n_contr, mesh] = getParams(testCase);
                
                [h, anaDiff] = testCase.h();
                
                
                opts_ = testCase.dyn.solver.opts;
                
                func = @() testCase.h();
                numDiff = testCase.numDiff_nD_AllT(func);
                testCase.assertSize(anaDiff, size(numDiff));
                testCase.assertLessThan(max(abs(anaDiff - numDiff)), opts_.RelTol);
            end
            
            %% ToDo: Eulerabfrage
            %         function testhDD(obj)
            %
            %             % TESTHDD This methods derives numerically obj.hd and compares
            %             % it with obj.hDD
            %             n_intervals = 3;
            %             obj.setupTest(n_intervals);
            %             [n_int, n_state, n_contr, mesh, n_var] = getParams(obj);
            %
            %             opts_ = obj.dyn.solver.opts;
            %             func = @() obj.gethD;
            %             numDiff = obj.numDiff_nxnD_AllT(func);
            %             anaDiff = obj.hDD();
            %
            %             size_nDiff_i = (n_var) * (n_intervals +1 );
            %             for i = 1:(n_int * n_state)
            %                 numDiff_i = reshape(numDiff(i,:,:), [size_nDiff_i size_nDiff_i]);
            %                 obj.assertSize(anaDiff{1}, size(numDiff_i));
            %                 obj.assertLessThan(max(abs(anaDiff{i} - numDiff_i)), opts_.RelTol * 9);
            %             end
            %         end
        end
    end
