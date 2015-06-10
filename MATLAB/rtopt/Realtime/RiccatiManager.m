classdef RiccatiManager <  TestEnv
    % RICCATIMANAGER Central class to perform a riccati solution
    
    properties(Access=private)
        robot;
        A;
        B;
        M;
        P;
        Q;
        R;
        
        D;
        
        nabla_s_star;
        nabla_lambda;
        nabla_q;
        
        nabla_mu;
        
        delta_lambda;
        delta_s;
        delta_q;
        
        delta_mu;
        
        horizon;
        n_lambda;
        n_state;
        n_contr;
        n_mu; % Das ist eine Cell, da n_mu die aktiven Constraints in jedem Zeitschritt unterschiedlich sind.
        n_var;
    end
    
    properties
        delta;
        delta_lambda;
        delta_s;
        delta_q;
        
        delta_mu;
    end
    
    
    methods
        function o = RiccatiManager(varargin)
            
            if(nargin == 0)
                global TEST;
                if ~(~isempty(TEST) && TEST == true)
                    error('wrong number of inputs');
                else
                    o.horizon = 20;
                    o.robot = Quadrocopter();
                end
            elseif(nargin == 1)
                o.horizon = varargin{1};
                o.robot = Quadrocopter();
            elseif(nargin == 2)
                o.horizon = varargin{1};
                if (isa(varargin{2},'Model'))
                    o.robot = varargin{2};
                end
            else
                error('wrong number of inputs');
            end
            
            o.initialize();
        end
        
        function doStep(o,i, LDD_i, LD_i, n_mu_i )
            %TODO: debug
            o.n_mu{i} = n_mu_i;
            o.Q{i} = LDD_i(1:o.n_state,1:o.n_state);
            o.nabla_lambda{i} = LD_i(1:o.n_lambda);
            o.n_var{i} = o.n_lambda + o.n_state + o.n_contr + n_mu_i;
            
            if(i == o.horizon+1)
                %Dann ist n_mu_i = 0 da dann keine controls vorhanden sind
                o.P{i} = o.Q{i};
                o.nabla_s_star{i} = LD_i(o.n_lambda+1:o.n_lambda + o.n_state);
                
            else
                o.M{i} = LDD_i(1:o.n_state,o.n_state +1 :o.n_state + o.n_contr);
                o.R{i} = LDD_i(o.n_state +1 :o.n_state + o.n_contr, o.n_state+1 :o.n_state + o.n_contr);
                o.A{i} = LDD_i(o.n_state + o.n_contr + n_mu_i + 1: o.n_var{i}, 1:o.n_state);
                o.B{i} = LDD_i(o.n_state + o.n_contr + n_mu_i + 1: o.n_var{i}, o.n_state +1 : o.n_state + o.n_contr);
                o.nabla_q{i} = LD_i(o.n_lambda + o.n_state +1  : o.n_lambda + o.n_state + o.n_contr);
                if(n_mu_i == 0)
                    % Benutze Regel aus Riccati_woConstr
                    
                    o.P{i} = o.Q{i} + (o.A{i}' * o.P{i+1} * o.A{i}) - ...
                        (o.M{i} + o.A{i}' * o.P{i+1} * o.B{i}) * ...
                        ((o.R{i} + o.B{i}' * o.P{i+1}  * o.B{i}) \ ...
                        (o.M{i}'  + o.B{i}' * o.P{i+1} * o.A{i}));
                    
                    o.nabla_s_star{i} = LD_i(o.n_lambda +1 : o.n_lambda + o.n_state) + ...
                        o.A{i}' * o.P{i+1} * o.nabla_lambda{i+1} + ...
                        o.A{i}' * o.nabla_s_star{i+1} - ...
                        (o.M{i} + o.A{i}' * o.P{i+1} * o.B{i}) * ...
                        (( o.R{i} +  o.B{i}' * o.P{i+1} *   o.B{i}) ...
                        \ ( o.nabla_q{i} + o.B{i}' * o.P{i+1} *  o.nabla_lambda{i+1} + o.B{i}' * o.nabla_s_star{i+1}));
                    %TODO: hier kann man o.B{i}' und auch o.A{i}' ausklammern, machen wir
                    %aber erst, wenn der Test geht.
                else
                    % Benutze neue Regel mit Constraints
                    o.D{i} = LDD_i(o.n_state + o.n_contr + 1 : o.n_state + o.n_contr + n_mu_i ,  o.n_state +1 : o.n_state + o.n_contr);
                    o.nabla_mu{i} = LD_i( o.n_lambda + o.n_state + o.n_contr +1 : o.n_lambda + o.n_state + o.n_contr + n_mu_i);
                    
                    z4 = [o.M{i}'  + o.B{i}' * o.P{i+1} * o.A{i}; zeros(n_mu_i, o.n_state)];
                    z3 = [ o.R{i} +  o.B{i}' * o.P{i+1} *   o.B{i}, o.D{i}' ;...
                        o.D{i} , zeros(n_mu_i)];
                    z2 = [o.M{i} + o.A{i}' * o.P{i+1} * o.B{i}, zeros(o.n_state, n_mu_i)];
                    z1 =  o.Q{i} + (o.A{i}' * o.P{i+1} * o.A{i});
                    
                    o.P{i} = z1 - z2 * (z3 \ z4);
                    
                    z5 = [ o.nabla_q{i} + o.B{i}' * ( o.P{i+1} * o.nabla_lambda{i+1} + o.nabla_s_star{i+1}) ; ...
                        o.nabla_mu{i}];
                    
                    z6 = LD_i(o.n_lambda +1 : o.n_lambda + o.n_state) + ...
                        o.A{i}' * (o.P{i+1} * o.nabla_lambda{i+1} + o.nabla_s_star{i+1});
                    
                    o.nabla_s_star{i} = z6 - z2 * (z3 \ z5);
                    
                end
            end
        end
        
        function solveStep(o,i)
            %TODO: implement
        end
        
        
        
    end
    
    methods(Test)
        
        
        
    end
    
    methods
        
        function initialize(o)
            %Initialize storage
            o.n_lambda = o.robot.n_state;
            o.n_state = o.robot.n_state;
            o.n_contr = o.robot.n_contr;
            
            
            o.A = cell(o.horizon,1);
            o.B = cell(o.horizon,1);
            o.M = cell(o.horizon,1);
            o.P = cell(o.horizon+1,1);
            o.Q = cell(o.horizon+1,1);
            o.R = cell(o.horizon,1);
            o.C = cell(o.horizon,1);
            o.D = cell(o.horizon,1);
            
            o.nabla_s_star = cell(o.horizon+1, 1);
            o.nabla_lambda = cell(o.horizon+1,1);
            o.nabla_q = cell(o.horizon ,1);
            o.nabla_mu = cell(o.horizon, 1);
            
            o.delta = zeros( o.n_var * (o.horizon +1 ) - (o.n_mu +  o.n_contr), 1);
            
            
            o.delta_lambda = cell(o.horizon +1, 1);
            o.delta_s = cell(o.horizon +1, 1);
            o.delta_q = cell(o.horizon, 1);
            o.delta_mu = cell(o.horizon,1);
        end
    end
end