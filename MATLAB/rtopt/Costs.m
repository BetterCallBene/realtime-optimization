classdef Costs < TestEnv
    % COST providing cost function based on control only
    %   Approximate the integral over the squared two-norm of the controls
    %   by a piecewise constant function using the discretized control values,
    %   i.e., the cost value is a sum of the the squared control values
    %   multiplied by the interval widths.
    
    properties
        dyn; % Dynamik
    end
    
    properties  (Dependent)
        vec;  % optimization vector combining all control and state values
    end
    
    methods
        %constructor
        function cC = Costs(varargin)
            % constructor based on one input value: a classOCPparam object
            if(nargin == 0)
                global TEST;                
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                end
                % constructor based on two input values
                % a classDyn element and a classOCPparam element
            elseif (nargin == 1)
                mc = metaclass(varargin{1});
                if strcmp(mc.SuperclassList(1).SuperclassList(1).Name, 'Dyn')
                    cC.dyn = varargin{1};
                else
                    error('wrong class type for problem parameter');
                end
            else
                error('wrong number of inputs');
            end
        end
        
        %set methods
        function set.vec(obj, vec_)
            % set new input vector
            obj.dyn.backdoor_vec = vec_;
        end
        
        function vec = get.vec(obj)
            vec = obj.dyn.vec;
        end
        
        % other functions
        function c_val = get_cost(obj)
            % compute the cost value
            control = obj.dyn.contr;
            mesh    = obj.dyn.environment.mesh;
            c_val   = .5*sum((control(:,1:end-1).^2)*mesh');
        end
        
        function cD_val = get_costD(obj)
            % compute the Jacobian of the cost value
            [n_int, n_tp, n_state, n_contr] = getParams(obj);
            
            control     = obj.dyn.contr;
            mesh        = obj.dyn.environment.mesh;
            
            rind        = zeros(1,(n_tp-1)*n_contr);
            vind        = zeros(1,(n_tp-1)*n_contr);
            for i = 1:n_tp-1
                rind((i-1)*n_contr+1:i*n_contr) = ...
                    (i-1)*(n_state+n_contr)+n_state + 1:...
                    i*(n_state+n_contr);
                vind((i-1)*n_contr+1:i*n_contr) = ...
                    mesh(i)*control(:,i)';
            end
            
            cD_val        = sparse(rind,ones(1,(n_tp-1)*n_contr),vind,...
                n_tp*(n_state+n_contr),1);
        end
        
        function cDD_val = get_costDD(obj)
            %compute the Hessian of the cost value
            [n_int, n_tp, n_state, n_contr] = getParams(obj);
            
            mesh        = obj.dyn.environment.mesh;
            
            rind        = zeros(1,(n_tp-1)*n_contr);
            vind        = zeros(1,(n_tp-1)*n_contr);
            for i = 1:n_tp-1
                rind((i-1)*n_contr+1:i*n_contr) = ...
                    (i-1)*(n_state+n_contr)+n_state + 1:...
                    i*(n_state+n_contr);
                vind((i-1)*n_contr+1:i*n_contr) = ...
                    mesh(i)*ones(1,n_contr);
            end
            
            cDD_val        = {sparse(rind,rind,vind,...
                n_tp*(n_state+n_contr),n_tp*(n_state+n_contr))};
        end
        
        function [n_int, n_tp, n_state, n_contr] = getParams(obj)
            n_int       = obj.dyn.environment.n_intervals;
            n_tp        = obj.dyn.environment.n_timepoints;
            n_state     = obj.dyn.robot.n_state;
            n_contr     = obj.dyn.robot.n_contr;
        end
        % general type get functions -> interace for testing
        function f = get_func(obj)
            f = get_cost(obj);
        end
        
        function g = get_jac(obj)
            g = get_costD(obj);
        end
        
        function H = get_hess(obj)
            H = get_costDD(obj);
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
            obj.vec = obj.dyn.vec;
            
        end
        
        function [vec_old, n, m, n_timepoints, dyn] = setup(obj, func)
            dyn = obj.dyn;
            vec_old = obj.vec;
            n_timepoints = obj.dyn.environment.n_timepoints;
            n = obj.dyn.robot.n_var;
            m = size(func());
            m = m(1);
        end
    end
    
    methods(Test)
        function testget_costD(obj)
            n_intervals = 50;
            obj.setupTest(n_intervals);
            
            func = @() obj.get_cost();
            anaDiff = obj.get_costD()';
            numDiff = obj.numDiff_nD_AllT(func);
            
            obj.assertSize(anaDiff, size(numDiff));
            obj.assertLessThan(max(abs(numDiff - anaDiff)),obj.tol);
            
        end
        
        function testget_costDD(obj)
            n_intervals = 50;
            obj.setupTest(n_intervals);
            
            func = @() obj.get_costD();
            anaDiff = obj.get_costDD();
            numDiff = obj.numDiff_nD_AllT(func);
            
            obj.assertSize(anaDiff, [1 1]);
            anaDiff = anaDiff{1};
            obj.assertSize(anaDiff, size(numDiff));
            obj.assertLessThan(max(abs(numDiff - anaDiff)), obj.tol);
        end
    end
end